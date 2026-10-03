"""
Wrong-Way River: the wrong-way raindrop scanner.  See README.md for setup, data, and usage.

Finds land cells whose water travels a long way to its endpoint compared with the
straight-line distance to the nearest ocean coast.

Run:  wrong-way-river     (after `pip install -e .`)   or   python wrong_way_river.py

All settings are environment variables (all optional):
  WRONG_WAY_TIF              flow-direction GeoTIFF (default hyd_glo_dir_30s.tif)
  WRONG_WAY_BBOX             "west,south,east,north": only read this part of the file.
                             Must contain whole basins AND the coasts near them.
                             Default: the whole file.
  WRONG_WAY_REGION           "west,south,east,north": the whole file is still used for the
                             flow paths and coasts, but only start cells inside this box are
                             ranked (e.g. an "Asia only" leaderboard). Default: everywhere.
  WRONG_WAY_UNITS            mi (default) or km, for everything that is printed or written
  WRONG_WAY_TOP_N            results per leaderboard (default 20)
  WRONG_WAY_PARETO_N         max entries per Pareto frontier list (default 100)
  WRONG_WAY_KML_TOP          results per non-Pareto leaderboard drawn in the KML
                             (default 10, 0 = no KML). Pareto entries are always all drawn.

Distance settings: give each as <NAME>_MI or <NAME>_KM (miles win if both are set).
  WRONG_WAY_DECLUSTER        keep one result per this radius (default 50 km = 31 mi)
  WRONG_WAY_PARETO_DECLUSTER same, for the Pareto lists (default 15 km = 9 mi)
  WRONG_WAY_BASIN_SEP        two results at the same river mouth are both listed only if
                             this far apart; 0 = no basin rule (default 1000 km = 621 mi)
  WRONG_WAY_MOUTH_MERGE      big river mouths closer than this count as one mouth, so a
                             delta's many channels are one basin (default 100 km = 62 mi)
  WRONG_WAY_BAND             near-coast band used to prune the search (default 100 km)
  WRONG_WAY_COAST_CUTOFF     a flow path ending within this distance of the coast counts as
                             ocean-draining, farther = endorheic (default 30 km = 19 mi)

Other:
  WRONG_WAY_FORCE_OCEAN      "lat,lon;lat,lon": treat the river mouth nearest each point as
                             ocean-draining (overrides the cutoff)
  WRONG_WAY_CACHE_DIR        cache folder (default wrong_way_cache)
  WRONG_WAY_REFRESH          set to anything to ignore all caches and recompute
"""
import os
import sys
import json
import time
import pickle
import hashlib
from xml.sax.saxutils import escape

import numpy as np
import pandas as pd
import rasterio
from affine import Affine
from rasterio.windows import from_bounds
from numba import njit
from scipy.ndimage import binary_dilation, label
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

__version__ = "0.5.0"

KM_PER_MI = 1.609344
MI_PER_KM = 1 / KM_PER_MI


def dist_env(base, default_km):
    """Read a distance setting given as <base>_MI or <base>_KM; returns kilometres."""
    mi = os.environ.get(base + "_MI")
    if mi:
        return float(mi) * KM_PER_MI
    km = os.environ.get(base + "_KM")
    return float(km) if km else default_km


def parse_box(text):
    if not text:
        return None
    box = tuple(float(v) for v in text.split(","))
    if len(box) != 4:
        raise SystemExit("boxes must be 'west,south,east,north'")
    return box


DIR_PATH = os.environ.get("WRONG_WAY_TIF", "hyd_glo_dir_30s.tif")
BBOX_ENV = os.environ.get("WRONG_WAY_BBOX")
REGION_ENV = os.environ.get("WRONG_WAY_REGION", "")
REGION = parse_box(REGION_ENV)
UNITS = os.environ.get("WRONG_WAY_UNITS", "mi").lower()
if UNITS not in ("mi", "km"):
    raise SystemExit("WRONG_WAY_UNITS must be 'mi' or 'km'")
UNIT_FACTOR = MI_PER_KM if UNITS == "mi" else 1.0
TOP_N = int(os.environ.get("WRONG_WAY_TOP_N", "20"))
PARETO_N = int(os.environ.get("WRONG_WAY_PARETO_N", "100"))
KML_TOP = int(os.environ.get("WRONG_WAY_KML_TOP", "10"))
DECLUSTER_KM = dist_env("WRONG_WAY_DECLUSTER", 50.0)
PARETO_DECLUSTER_KM = dist_env("WRONG_WAY_PARETO_DECLUSTER", 15.0)
BASIN_SEP_KM = dist_env("WRONG_WAY_BASIN_SEP", 1000.0)
MOUTH_MERGE_KM = dist_env("WRONG_WAY_MOUTH_MERGE", 100.0)
BAND_KM = dist_env("WRONG_WAY_BAND", 100.0)
COAST_CUTOFF_KM = dist_env("WRONG_WAY_COAST_CUTOFF", 30.0)
FORCE_OCEAN = os.environ.get("WRONG_WAY_FORCE_OCEAN", "")
CACHE_DIR = os.environ.get("WRONG_WAY_CACHE_DIR", "wrong_way_cache")
REFRESH = bool(os.environ.get("WRONG_WAY_REFRESH"))

SEP = BASIN_SEP_KM if BASIN_SEP_KM > 0 else None
MIN_STRAIGHT_KM = 5.0       # floor for the ratio
MIN_DETOUR_KM = 100.0       # ignore cells whose detour can't exceed this
OCEAN_MIN_CELLS = 5000      # nodata blobs bigger than this are "ocean" (used for the coast)
MOUTH_MIN_CELLS = 20_000    # only river mouths draining at least this many cells get merged
BORDER_MAX_KM = 150.0       # upper end of the "borderline mouths" table
FORCE_RADIUS_KM = 25.0
CHUNK = 500_000
R = 6371.0

# leaderboard titles and short tags, in priority order (most salient first)
T_OPAR = "Ocean-draining: Pareto frontier (approximate)"
T_ODET = "Ocean-draining: biggest absolute detour"
T_ORAT = "Ocean-draining: biggest ratio"
T_EPAR = "Endorheic: Pareto frontier (approximate)"
T_EDET = "Endorheic: biggest absolute detour"
T_ERAT = "Endorheic: biggest ratio"
PARETO_TITLES = (T_OPAR, T_EPAR)
LIST_TAGS = [(T_OPAR, "Pareto", 0), (T_ODET, "Detour", 1), (T_ORAT, "Ratio", 1),
             (T_EPAR, "Endo Pareto", 2), (T_EDET, "Endo Detour", 2), (T_ERAT, "Endo Ratio", 2)]
TIER_NAMES = {0: "Pareto frontier (ocean-draining)",
              1: "Other top results (ocean-draining)",
              2: "Endorheic (closed basins)"}

DIRS = {1: (0, 1), 2: (1, 1), 4: (1, 0), 8: (1, -1),
        16: (0, -1), 32: (-1, -1), 64: (-1, 0), 128: (-1, 1)}
DR = np.zeros(256, np.int64)
DC = np.zeros(256, np.int64)
VALID = np.zeros(256, np.bool_)
for k, (dr, dc) in DIRS.items():
    DR[k], DC[k], VALID[k] = dr, dc, True

MOUTH_GROUPS = {"u": None, "g": None}
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:6.0f}s] {msg}", flush=True)


def fmt_d(km):
    """Distance in the display units, e.g. '62 mi'."""
    return f"{km * UNIT_FACTOR:,.0f} {UNITS}"


def maps_link(lat, lon):
    return f"https://www.google.com/maps?q={lat:.5f},{lon:.5f}"


# --------------------------------------------------------------------------- flow walk
@njit(cache=True)
def flow_distance(fdir, land, dy_km, dx_km, DR, DC, VALID):
    """
    dist: km along the flow path to the terminal cell (-1 for non-land).
    term: flat index of the terminal cell the water ends in,
          -2 if the path leaves the box / is invalid, -1 for non-land.
    """
    H, W = fdir.shape
    dist = np.full((H, W), -1.0, np.float32)
    term = np.full((H, W), -1, np.int32)
    stack = np.empty(5_000_000, np.int32)
    cap = stack.shape[0]
    for sr in range(H):
        for sc in range(W):
            if not land[sr, sc] or dist[sr, sc] >= 0:
                continue
            n = 0
            r, c = sr, sc
            while True:
                if dist[r, c] >= 0:
                    break                                   # joined a solved path
                code = fdir[r, c]
                if code <= 0:                               # sink cell (outlet or inland)
                    dist[r, c] = 0.0
                    term[r, c] = r * W + c
                    break
                if code > 255 or not VALID[code]:
                    dist[r, c] = 0.0; term[r, c] = -2; break
                nr, nc = r + DR[code], c + DC[code]
                if nr < 0 or nr >= H or nc < 0 or nc >= W:
                    dist[r, c] = 0.0; term[r, c] = -2; break    # leaves the box
                if not land[nr, nc]:                        # steps into nodata
                    dist[r, c] = 0.0
                    term[r, c] = r * W + c
                    break
                if n >= cap:
                    dist[r, c] = 0.0; term[r, c] = -2; break    # absurdly long path
                stack[n] = r * W + c
                n += 1
                r, c = nr, nc
            pr, pc = r, c
            for i in range(n - 1, -1, -1):
                cr, cc = stack[i] // W, stack[i] % W
                dr, dc = abs(cr - pr), abs(cc - pc)
                if dr and dc:
                    step = np.hypot(dy_km, dx_km[cr])
                elif dr:
                    step = dy_km
                else:
                    step = dx_km[cr]
                dist[cr, cc] = dist[pr, pc] + step
                term[cr, cc] = term[pr, pc]
                pr, pc = cr, cc
    return dist, term


# --------------------------------------------------------------------------- geometry
def to_xyz(lat, lon):
    la = np.radians(np.asarray(lat, dtype=np.float64))
    lo = np.radians(np.asarray(lon, dtype=np.float64))
    return np.column_stack([R * np.cos(la) * np.cos(lo),
                            R * np.cos(la) * np.sin(lo),
                            R * np.sin(la)])


def km_to_chord(km):
    return 2 * R * np.sin(km / (2 * R))


def chord_to_km(d):
    d = np.asarray(d, dtype=np.float64)
    out = np.full(d.shape, np.inf)
    ok = np.isfinite(d)
    out[ok] = 2 * R * np.arcsin(np.clip(d[ok] / (2 * R), 0, 1))
    return out


def straight_for(tree, lats, lons, W, flat_idx, bound_km=None):
    """Great-circle km to the nearest coast cell; inf where beyond bound_km."""
    out = np.empty(len(flat_idx), np.float32)
    ub = np.inf if bound_km is None else km_to_chord(bound_km)
    for i in range(0, len(flat_idx), CHUNK):
        f = flat_idx[i:i + CHUNK]
        r, c = f // W, f % W
        d, _ = tree.query(to_xyz(lats[r], lons[c]), distance_upper_bound=ub, workers=-1)
        out[i:i + CHUNK] = chord_to_km(d)
    return out


# --------------------------------------------------------------------------- river-mouth groups
def build_mouth_groups(u, counts, lats, lons, W):
    """
    Group id per river mouth. Big mouths (drain >= MOUTH_MIN_CELLS cells) closer than
    MOUTH_MERGE_KM are merged, so a delta's many channels count as one basin.
    """
    group = u.copy()
    big = np.flatnonzero(counts >= MOUTH_MIN_CELLS)
    if MOUTH_MERGE_KM > 0 and len(big) > 1:
        tree = cKDTree(to_xyz(lats[u[big] // W], lons[u[big] % W]))
        pairs = tree.query_pairs(km_to_chord(MOUTH_MERGE_KM), output_type="ndarray")
        if len(pairs):
            n = len(big)
            m = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
            _, labels = connected_components(m, directed=False)
            t = pd.DataFrame({"label": labels, "i": big, "c": counts[big]})
            rep = t.loc[t.groupby("label").c.idxmax()].set_index("label").i
            group[big] = u[rep.loc[labels].values]
    return group


def group_of(ids):
    ids = np.asarray(ids)
    if MOUTH_GROUPS["u"] is None:
        return ids
    return MOUTH_GROUPS["g"][np.searchsorted(MOUTH_GROUPS["u"], ids)]


# --------------------------------------------------------------------------- grouping results
def decluster(df, key, radius_km=None, n=None, pool=500_000, basin_sep_km=None):
    """
    Greedy: take the best row by `key`, drop everything within radius_km of it, repeat.
    If basin_sep_km is set, a row is also skipped when an already-kept row is at the same
    river mouth (same mouth group) and closer than basin_sep_km.
    """
    radius_km = DECLUSTER_KM if radius_km is None else radius_km
    n = TOP_N if n is None else n
    d = df.nlargest(pool, key).reset_index(drop=True)
    if d.empty:
        return d
    xyz = to_xyz(d.lat.values, d.lon.values)
    tree = cKDTree(xyz)
    chord = km_to_chord(radius_km)
    sep_chord = km_to_chord(basin_sep_km) if basin_sep_km else None
    mids = group_of(d["mouth_id"].values) if sep_chord is not None else None
    removed = np.zeros(len(d), bool)
    by_mouth = {}
    keep = []
    for i in range(len(d)):
        if removed[i]:
            continue
        if sep_chord is not None:
            prev = by_mouth.get(int(mids[i]))
            if prev is not None and np.min(np.linalg.norm(xyz[prev] - xyz[i], axis=1)) < sep_chord:
                continue
            by_mouth.setdefault(int(mids[i]), []).append(i)
        keep.append(i)
        if len(keep) >= n:
            break
        removed[tree.query_ball_point(xyz[i], r=chord)] = True
    return d.loc[keep]


def nth_detour(d):
    """Detour of the N-th declustered result (0 if there are fewer than N)."""
    top = decluster(d, "detour_km", basin_sep_km=SEP)
    return float(top.detour_km.iloc[-1]) if len(top) >= TOP_N else 0.0


def build_stories(lists, top=None):
    """
    Merge entries that describe the same place across all leaderboards. Two entries are
    the same story if they are within the decluster radius of each other, or share a river
    mouth group and are within the basin separation. The first (most salient) list an
    entry appears in provides the representative cell.
    """
    radius_chord = km_to_chord(DECLUSTER_KM)
    sep_chord = km_to_chord(BASIN_SEP_KM) if SEP else None
    stories = []
    for title, tag, tier in LIST_TAGS:
        d = lists.get(title)
        if d is None or not len(d):
            continue
        rows = d if (title in PARETO_TITLES or top is None) else d.head(top)
        for rank, row in enumerate(rows.itertuples(), 1):
            x = to_xyz([row.lat], [row.lon])[0]
            gid = int(group_of(np.array([row.mouth_id]))[0])
            hit = None
            for s in stories:
                for mx, mg in s["pts"]:
                    dist = np.linalg.norm(x - mx)
                    if dist < radius_chord or (sep_chord is not None and mg == gid and dist < sep_chord):
                        hit = s
                        break
                if hit:
                    break
            if hit is None:
                hit = {"rep": row, "tier": tier, "ranks": [], "pts": [], "members": []}
                stories.append(hit)
            hit["ranks"].append(f"{tag} #{rank}")
            hit["pts"].append((x, gid))
            hit["members"].append((f"{tag} #{rank}", row))
    return stories


# --------------------------------------------------------------------------- caching steps
def compute_flow(ds, bbox, paths):
    """Read the raster, find the coast, walk every flow path, and save it all to the cache."""
    win = from_bounds(*bbox, transform=ds.transform)
    fdir = ds.read(1, window=win).astype(np.int16)
    tf = ds.window_transform(win)
    nodata = ds.nodata
    H, W = fdir.shape
    log(f"window: {H} x {W} cells")
    land = (fdir != nodata) if nodata is not None else (fdir >= 0)

    log("finding the ocean (for the coastline)...")
    lab, _ = label(~land)
    sizes = np.bincount(lab.ravel())
    big = sizes > OCEAN_MIN_CELLS
    big[0] = False
    ocean_mask = big[lab]
    del lab
    log(f"  ocean blobs: {int(big.sum())}")
    ocean_dil = binary_dilation(ocean_mask, structure=np.ones((3, 3), bool))
    cr, cc = np.nonzero(land & ocean_dil)
    del ocean_mask, ocean_dil

    lats = tf.f + (np.arange(H) + 0.5) * tf.e
    dy_km = np.radians(tf.a) * R
    dx_km = (np.radians(tf.a) * R * np.cos(np.radians(lats))).astype(np.float64)

    log("flow distances (the slow step; cached afterwards)...")
    fd, term = flow_distance(fdir, land, dy_km, dx_km, DR, DC, VALID)
    del fdir, land

    os.makedirs(os.path.dirname(paths["fd"]), exist_ok=True)
    np.save(paths["fd"], fd)
    np.save(paths["term"], term)
    np.savez(paths["coast"], cr=cr.astype(np.int32), cc=cc.astype(np.int32))
    json.dump({"H": H, "W": W, "transform": list(tf)[:6]}, open(paths["meta"], "w"))
    log("cached flow results")


def build_terminal_table(term, tree, lats, lons, W, ds, bbox):
    """One row per river mouth / sink: upstream cell count, distance to coast, HydroSHEDS code."""
    log("building the table of river mouths and sinks...")
    t = np.asarray(term)
    u, counts = np.unique(t[t >= 0], return_counts=True)
    n_trunc = int((t == -2).sum())
    tr, tc = u // W, u % W
    chord, _ = tree.query(to_xyz(lats[tr], lons[tc]), workers=-1)
    t_km = chord_to_km(chord)
    win = from_bounds(*bbox, transform=ds.transform)
    fdir = ds.read(1, window=win)
    tcode = fdir.ravel()[u].astype(np.int16)
    return u, counts, t_km, tcode, n_trunc


def print_terminals(title, idx, u, counts, t_km, tcode, lats, lons, W, limit):
    print(f"\n{title}")
    if not len(idx):
        print("  (none)")
        return
    print(f"     lat       lon       cells  {UNITS}_from_coast  HydroSHEDS_flag  map")
    for i in idx[np.argsort(counts[idx])[::-1][:limit]]:
        la, lo = lats[u[i] // W], lons[u[i] % W]
        flag = "outlet(0)" if tcode[i] == 0 else "inland-sink"
        print(f"  {la:8.3f} {lo:9.3f} {counts[i]:>11,} {t_km[i] * UNIT_FACTOR:>14.1f}  {flag:<15}  "
              f"{maps_link(la, lo)}")


def parse_force(text):
    pts = []
    for part in text.split(";"):
        part = part.strip()
        if part:
            la, lo = part.split(",")
            pts.append((float(la), float(lo)))
    return pts


# --------------------------------------------------------------------------- scoring
def score_cells(fd, term, u, t_ocean, tree, lats, lons, W):
    """
    Returns (df, T_used, unc_max).
    Phase 1: straight distance for cells within BAND_KM of the coast.
    Phase 2: inland cells only if their flow distance exceeds the current N-th best detour
             (detour = flow - straight <= flow, so smaller-flow cells cannot make the list).
    """
    valid = np.asarray(term) >= 0
    flat = np.flatnonzero(valid.ravel())
    del valid
    fl = np.asarray(fd).ravel()[flat]
    keep = fl > MIN_DETOUR_KM
    flat, fl = flat[keep], fl[keep]
    if REGION:
        w_, s_, e_, n_ = REGION
        rr, qq = flat // W, flat % W
        m = (lats[rr] >= s_) & (lats[rr] <= n_) & (lons[qq] >= w_) & (lons[qq] <= e_)
        flat, fl = flat[m], fl[m]
        log(f"region filter {REGION_ENV}: {len(flat):,} cells inside")
    tflat = np.asarray(term).ravel()[flat]
    cls_endo = ~t_ocean[np.searchsorted(u, tflat)]
    log(f"candidate cells (flow > {fmt_d(MIN_DETOUR_KM)}): {len(flat):,}")

    def build_df(idx):
        f = flat[idx]
        r, c = f // W, f % W
        d = pd.DataFrame({"cell_id": f.astype(np.int32),
                          "mouth_id": tflat[idx],
                          "lat": lats[r].astype(np.float32),
                          "lon": lons[c].astype(np.float32),
                          "flow_km": fl[idx],
                          "straight_km": straight[idx],
                          "endorheic": cls_endo[idx]})
        d["detour_km"] = d.flow_km - d.straight_km
        d["ratio"] = d.flow_km / d.straight_km.clip(lower=MIN_STRAIGHT_KM)
        return d[d.detour_km > MIN_DETOUR_KM]

    log(f"phase 1: straight distance for cells within {fmt_d(BAND_KM)} of the coast...")
    straight = straight_for(tree, lats, lons, W, flat, BAND_KM)
    done = np.isfinite(straight)
    log(f"  {int(done.sum()):,} cells in the band")
    df1 = build_df(np.flatnonzero(done))
    T = {False: nth_detour(df1[~df1.endorheic]), True: nth_detour(df1[df1.endorheic])}
    log(f"  current N-th best detour: ocean {fmt_d(T[False])}, endorheic {fmt_d(T[True])}")
    del df1

    thr = np.where(cls_endo, T[True], T[False]).astype(np.float32)
    need = np.flatnonzero(~done & (fl > thr))
    log(f"phase 2: {len(need):,} inland cells could still beat that; computing them...")
    straight[need] = straight_for(tree, lats, lons, W, flat[need], None)

    computed = np.isfinite(straight)
    log(f"computed {int(computed.sum()):,} of {len(flat):,} candidate cells")
    df = build_df(np.flatnonzero(computed))
    unc = ~computed
    unc_max = {c: (float(fl[unc & (cls_endo == c)].max()) if (unc & (cls_endo == c)).any() else 0.0)
               for c in (False, True)}
    return df, T, unc_max


# --------------------------------------------------------------------------- output
def add_geo(d, g):
    """Add river-mouth and nearest-coast coordinates to a (small) results table."""
    d = d.copy()
    if not len(d):
        return d
    mid = d.mouth_id.values.astype(np.int64)
    d["mouth_lat"] = g["lats"][mid // g["W"]]
    d["mouth_lon"] = g["lons"][mid % g["W"]]
    _, ci = g["tree"].query(to_xyz(d.lat.values, d.lon.values), workers=-1)
    d["coast_lat"] = g["lats"][g["cr"][ci]]
    d["coast_lon"] = g["lons"][g["cc"][ci]]
    return d


def with_units(d):
    """Replace the internal *_km columns with columns in the display units."""
    names = ("flow", "straight", "detour")
    vals = {n: d[f"{n}_km"] * UNIT_FACTOR for n in names}
    d = d.drop(columns=[f"{n}_km" for n in names])
    for n in names:
        d[f"{n}_{UNITS}"] = vals[n]
    return d


def show(title, d, lead=None):
    print(f"\n=== {title} ===")
    if not len(d):
        print("(none)")
        return
    d = with_units(d)
    d["map"] = [maps_link(a, b) for a, b in zip(d.lat, d.lon)]
    d["mouth"] = [maps_link(a, b) for a, b in zip(d.mouth_lat, d.mouth_lon)]
    cols = ([lead] if lead else []) + ["lat", "lon", f"flow_{UNITS}", f"straight_{UNITS}",
                                      f"detour_{UNITS}", "ratio", "endorheic", "map", "mouth"]
    print(d[cols].round(2).to_string(index=False))


def trace_path(fd, r, c, dy_km, dx_km, H, W, max_steps=300_000):
    """
    Follow the water downstream using the cached flow-distance grid: the next cell is the
    neighbor whose distance plus the step to it equals this cell's distance.
    """
    pts = [(r, c)]
    cur = float(fd[r, c])
    for _ in range(max_steps):
        if cur <= 0:
            break
        best = None
        for ddr in (-1, 0, 1):
            for ddc in (-1, 0, 1):
                if ddr == 0 and ddc == 0:
                    continue
                nr, nc = r + ddr, c + ddc
                if nr < 0 or nr >= H or nc < 0 or nc >= W:
                    continue
                v = float(fd[nr, nc])
                if v < 0 or v >= cur:
                    continue
                if ddr and ddc:
                    step = float(np.hypot(dy_km, dx_km[r]))
                elif ddr:
                    step = dy_km
                else:
                    step = float(dx_km[r])
                err = abs(v + step - cur)
                if best is None or err < best[0]:
                    best = (err, nr, nc, v)
        if best is None or best[0] > 0.05:
            break
        _, r, c, cur = best
        pts.append((r, c))
    return pts


# KML colors are aabbggrr. Pareto = bold red path + green straight line; others are muted.
KML_STYLES = {
    0: ("ff0000ff", "ff00c800", 4),   # Pareto frontier, ocean-draining
    1: ("cc3399ff", "ccffc060", 2),   # other ocean-draining results
    2: ("aa808080", "aac0c0c0", 2),   # endorheic
}


def write_kml(path, stories, fd, g):
    H, W, lats, lons = g["H"], g["W"], g["lats"], g["lons"]
    out = ['<?xml version="1.0" encoding="UTF-8"?>',
           '<kml xmlns="http://www.opengis.net/kml/2.2"><Document>',
           '<name>Wrong-way raindrops</name>']
    for tier, (fc, sc, width) in KML_STYLES.items():
        out.append(f'<Style id="f{tier}"><LineStyle><color>{fc}</color><width>{width}</width></LineStyle></Style>')
        out.append(f'<Style id="s{tier}"><LineStyle><color>{sc}</color><width>{width}</width></LineStyle></Style>')
    n_bad = 0
    for tier in (0, 1, 2):
        group = [s for s in stories if s["tier"] == tier]
        if not group:
            continue
        out.append(f"<Folder><name>{escape(TIER_NAMES[tier])}</name>")
        for s in group:
            row = s["rep"]
            r, c = divmod(int(row.cell_id), W)
            mr, mc = divmod(int(row.mouth_id), W)
            pts = trace_path(fd, r, c, g["dy_km"], g["dx_km"], H, W)
            complete = pts[-1] == (mr, mc)
            n_bad += (not complete)
            step = max(1, len(pts) // 800)
            sel = pts[::step]
            if sel[-1] != pts[-1]:
                sel.append(pts[-1])
            flow_xy = " ".join(f"{lons[pc]:.5f},{lats[pr]:.5f},0" for pr, pc in sel)
            ranks = ", ".join(s["ranks"])
            extra = f" (+{len(s['ranks']) - 1})" if len(s["ranks"]) > 1 else ""
            short = (f"{s['ranks'][0]}{extra} \u00b7 {row.detour_km * UNIT_FACTOR:,.0f} {UNITS}"
                     f" \u00b7 {row.ratio:,.0f}x")
            title = f"{ranks}: {row.detour_km * UNIT_FACTOR:,.0f} {UNITS} detour, {row.ratio:,.0f}x"
            members = "<br>".join(
                f"{lab}: {m.detour_km * UNIT_FACTOR:,.0f} {UNITS} detour, {m.ratio:,.1f}x "
                f"(<a href=\"{maps_link(m.lat, m.lon)}\">location</a>)" for lab, m in s["members"])
            where = "ends in an inland sink" if row.endorheic else "reaches the ocean"
            desc = (f"<![CDATA[Appears in: {members}<br><br>"
                    f"Representative path: {row.flow_km * UNIT_FACTOR:,.0f} {UNITS} by water, "
                    f"{row.straight_km * UNIT_FACTOR:,.0f} {UNITS} straight to the nearest coast. "
                    f"Water {where}.{'' if complete else ' (path trace stopped early)'}<br>"
                    f"<a href=\"{maps_link(row.mouth_lat, row.mouth_lon)}\">Mouth in Google Maps</a>]]>")
            out.append(f"<Folder><name>{escape(title)}</name>")
            out.append(f"<Placemark><name>{escape(short)}</name><description>{desc}</description>"
                       f"<Point><coordinates>{row.lon:.5f},{row.lat:.5f},0</coordinates></Point></Placemark>")
            out.append(f"<Placemark><name>Flow path to the mouth</name><styleUrl>#f{tier}</styleUrl>"
                       f"<LineString><tessellate>1</tessellate><altitudeMode>clampToGround</altitudeMode>"
                       f"<coordinates>{flow_xy}</coordinates></LineString></Placemark>")
            out.append(f"<Placemark><name>Straight line to the nearest coast</name><styleUrl>#s{tier}</styleUrl>"
                       f"<LineString><tessellate>1</tessellate><altitudeMode>clampToGround</altitudeMode>"
                       f"<coordinates>{row.lon:.5f},{row.lat:.5f},0 {row.coast_lon:.5f},{row.coast_lat:.5f},0"
                       f"</coordinates></LineString></Placemark>")
            out.append(f"<Placemark><name>River mouth / sink</name>"
                       f"<Point><coordinates>{row.mouth_lon:.5f},{row.mouth_lat:.5f},0</coordinates></Point></Placemark>")
            out.append("</Folder>")
        out.append("</Folder>")
    out.append("</Document></kml>")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(out))
    log(f"wrote {path}: {len(stories)} places" + (f" ({n_bad} path traces stopped early)" if n_bad else ""))


# --------------------------------------------------------------------------- main
def require_input_file():
    """Exit with a friendly message if the flow-direction raster is missing or unreadable."""
    if not os.path.isfile(DIR_PATH):
        print("\nCan't find the flow-direction file.\n\n"
              f"  Looked for:   {os.path.abspath(DIR_PATH)}\n"
              f"  Running from: {os.getcwd()}\n\n"
              "This script needs the HydroSHEDS flow-direction raster (30 arc-second, global\n"
              "GeoTIFF), which is too large to ship with the code. To fix this:\n"
              "  1. Download it. The README's \"Data\" section explains which file to get and where\n"
              "     (start at https://www.hydrosheds.org/products/hydrosheds).\n"
              "  2. Put it in the folder you run from, named hyd_glo_dir_30s.tif, or point to it:\n"
              "       PowerShell:  $env:WRONG_WAY_TIF = \"C:\\path\\to\\hyd_glo_dir_30s.tif\"\n"
              "       cmd:         set WRONG_WAY_TIF=C:\\path\\to\\hyd_glo_dir_30s.tif\n"
              "       bash:        WRONG_WAY_TIF=/path/to/hyd_glo_dir_30s.tif wrong-way-river\n",
              file=sys.stderr)
        sys.exit(1)
    try:
        with rasterio.open(DIR_PATH):
            pass
    except rasterio.errors.RasterioIOError as err:
        print(f"\nFound {os.path.abspath(DIR_PATH)} but couldn't read it as a raster ({err}).\n"
              "It should be the HydroSHEDS flow-direction GeoTIFF (for example hyd_glo_dir_30s.tif),\n"
              "unzipped. See the \"Data\" section of the README.\n", file=sys.stderr)
        sys.exit(1)


def main():
    if any(a in ("-h", "--help") for a in sys.argv[1:]):
        print(__doc__)
        return
    if "--version" in sys.argv[1:]:
        print(f"wrong-way-river {__version__}")
        return
    require_input_file()
    with rasterio.open(DIR_PATH) as ds:
        bbox = parse_box(BBOX_ENV) or tuple(ds.bounds)
        log(f"file: {DIR_PATH}\n         bbox (w,s,e,n): {bbox}")

        st = os.stat(DIR_PATH)
        # same key format as earlier versions, so an existing flow cache is reused
        key = hashlib.md5(
            f"{os.path.abspath(DIR_PATH)}|{st.st_size}|{st.st_mtime}|{bbox}|{OCEAN_MIN_CELLS}|v3"
            .encode()).hexdigest()[:12]
        cdir = os.path.join(CACHE_DIR, key)
        paths = {"fd": os.path.join(cdir, "fd.npy"), "term": os.path.join(cdir, "term.npy"),
                 "coast": os.path.join(cdir, "coast.npz"), "meta": os.path.join(cdir, "meta.json"),
                 "terminals": os.path.join(cdir, "terminals.npz")}
        flow_ok = all(os.path.exists(paths[k]) for k in ("fd", "term", "coast", "meta"))

        if REFRESH or not flow_ok:
            compute_flow(ds, bbox, paths)
        else:
            log(f"using cached flow results in {cdir}")

        meta = json.load(open(paths["meta"]))
        H, W = meta["H"], meta["W"]
        tf = Affine(*meta["transform"])
        lats = tf.f + (np.arange(H) + 0.5) * tf.e
        lons = tf.c + (np.arange(W) + 0.5) * tf.a
        fd = np.load(paths["fd"], mmap_mode="r")
        term = np.load(paths["term"], mmap_mode="r")
        z = np.load(paths["coast"])
        cr, cc = z["cr"], z["cc"]
        log(f"coast cells: {len(cr):,}; building KD-tree...")
        tree = cKDTree(to_xyz(lats[cr], lons[cc]))

        if os.path.exists(paths["terminals"]) and not REFRESH:
            z = np.load(paths["terminals"])
            u, counts, t_km, tcode, n_trunc = z["u"], z["counts"], z["t_km"], z["tcode"], int(z["n_trunc"])
        else:
            u, counts, t_km, tcode, n_trunc = build_terminal_table(term, tree, lats, lons, W, ds, bbox)
            np.savez(paths["terminals"], u=u, counts=counts, t_km=t_km, tcode=tcode, n_trunc=n_trunc)

    # ---- classify river mouths: ocean vs endorheic ---------------------------------
    log(f"land cells with a usable endpoint: {int(counts.sum()):,}; truncated: {n_trunc:,}")
    t_ocean = t_km <= COAST_CUTOFF_KM
    forced = parse_force(FORCE_OCEAN)
    if forced:
        ttree = cKDTree(to_xyz(lats[u // W], lons[u % W]))
        for la, lo in forced:
            d, i = ttree.query(to_xyz([la], [lo])[0])
            if chord_to_km(d) <= FORCE_RADIUS_KM:
                t_ocean[i] = True
                log(f"forced ocean-draining: mouth at {lats[u[i] // W]:.3f}, {lons[u[i] % W]:.3f} "
                    f"({int(counts[i]):,} cells)")
            else:
                log(f"WARNING: no river mouth within {fmt_d(FORCE_RADIUS_KM)} of {la}, {lo}; ignored")
    log(f"river mouths/sinks: {len(u):,} ({int(t_ocean.sum()):,} ocean, {int((~t_ocean).sum()):,} endorheic); "
        f"cells draining to ocean {int(counts[t_ocean].sum()):,}, to sinks {int(counts[~t_ocean].sum()):,}")
    group = build_mouth_groups(u, counts, lats, lons, W)
    MOUTH_GROUPS["u"], MOUTH_GROUPS["g"] = u, group
    log(f"merged nearby big river mouths: {len(np.unique(group)):,} basin groups from {len(u):,} mouths")
    endo_i = np.flatnonzero(~t_ocean)
    print_terminals("Biggest endorheic sinks by upstream cell count (sanity check):",
                    endo_i, u, counts, t_km, tcode, lats, lons, W, 10)
    border = np.flatnonzero((~t_ocean) & (t_km <= BORDER_MAX_KM))
    print_terminals(f"Borderline: classified endorheic but within {fmt_d(BORDER_MAX_KM)} of the coast "
                    f"(candidates for WRONG_WAY_FORCE_OCEAN or a larger cutoff):",
                    border, u, counts, t_km, tcode, lats, lons, W, 15)

    # ---- scored cells (cached) ------------------------------------------------------
    sk = hashlib.md5(f"{key}|{COAST_CUTOFF_KM}|{FORCE_OCEAN}|{BAND_KM}|{MIN_DETOUR_KM}|"
                     f"{MIN_STRAIGHT_KM}|{REGION_ENV}".encode()).hexdigest()[:10]
    spath = os.path.join(cdir, f"scored_{sk}.pkl")
    scored = None
    if os.path.exists(spath) and not REFRESH:
        log("loading cached scored cells...")
        with open(spath, "rb") as fh:
            scored = pickle.load(fh)
        df = scored["df"]
        T_new = {c: nth_detour(df[df.endorheic == c]) for c in (False, True)}
        if all(T_new[c] >= scored["T"][c] - 1e-3 for c in (False, True)):
            log("  cache is sufficient for the current grouping settings")
        else:
            log("  current grouping settings need cells the cache skipped; rescoring")
            scored = None
    if scored is None:
        df, T_used, unc_max = score_cells(fd, term, u, t_ocean, tree, lats, lons, W)
        scored = {"df": df, "T": T_used, "unc_max": unc_max}
        with open(spath, "wb") as fh:
            pickle.dump(scored, fh, protocol=pickle.HIGHEST_PROTOCOL)
        log("cached scored cells")
    df, unc_max = scored["df"], scored["unc_max"]

    # ---- leaderboards -----------------------------------------------------------------
    g = {"H": H, "W": W, "lats": lats, "lons": lons, "tree": tree, "cr": cr, "cc": cc,
         "dy_km": float(np.radians(tf.a) * R),
         "dx_km": (np.radians(tf.a) * R * np.cos(np.radians(lats))).astype(np.float64)}

    def dec(d, key_, n=None):
        return decluster(d, key_, n=n, basin_sep_km=SEP)

    def pareto(d, title):
        """Cells no other cell beats on both detour and ratio, thinned to distinct places."""
        s_ = d.sort_values("detour_km", ascending=False)
        front = s_[s_.ratio.cummax().diff().fillna(1) > 0]
        top = decluster(front, "detour_km", radius_km=PARETO_DECLUSTER_KM, n=PARETO_N + 1)
        if len(top) > PARETO_N:
            print(f"\nNOTE: {title} has more than {PARETO_N} entries; keeping the {PARETO_N} with the "
                  f"biggest detour. Raise WRONG_WAY_PARETO_N to see more.")
            top = top.head(PARETO_N)
        return top

    ocean_df, endo_df = df[~df.endorheic], df[df.endorheic]
    lists = {T_ODET: dec(ocean_df, "detour_km"), T_ORAT: dec(ocean_df, "ratio"),
             T_OPAR: pareto(ocean_df, T_OPAR),
             T_EDET: dec(endo_df, "detour_km"), T_ERAT: dec(endo_df, "ratio"),
             T_EPAR: pareto(endo_df, T_EPAR)}
    lists = {k: add_geo(v, g) for k, v in lists.items()}
    for title in (T_ODET, T_ORAT, T_OPAR, T_EDET, T_ERAT, T_EPAR):
        show(title, lists[title])

    # ---- exactness report -------------------------------------------------------------
    print("\n=== Exactness check ===")
    print("Absolute-detour lists: exact (skipped cells had flow below the N-th best detour).")
    suffix = UNITS.upper()
    for endo, name, top in ((False, "ocean-draining", lists[T_ORAT]), (True, "endorheic", lists[T_ERAT])):
        bound = unc_max[endo] / BAND_KM
        rn = float(top.ratio.iloc[-1]) if len(top) >= TOP_N else 0.0
        if rn <= 0 and bound > 0:
            print(f"Ratio list, {name}: fewer than {TOP_N} results, and skipped cells farther than "
                  f"{fmt_d(BAND_KM)} from the coast could reach ratio {bound:.1f}. "
                  f"Raise WRONG_WAY_BAND_{suffix}.")
        elif rn >= bound:
            print(f"Ratio list, {name}: exact (N-th ratio {rn:.1f} >= best any skipped cell "
                  f"could have, {bound:.1f}).")
        else:
            need_km = unc_max[endo] / max(rn, 1e-9)
            print(f"Ratio list, {name}: may miss cells farther than {fmt_d(BAND_KM)} from the coast "
                  f"(skipped cells could reach ratio {bound:.1f} vs N-th ratio {rn:.1f}). "
                  f"Set WRONG_WAY_BAND_{suffix}={need_km * UNIT_FACTOR:.0f} to be sure.")

    # ---- consolidated highlights ------------------------------------------------------
    full = build_stories(lists)
    rows = []
    for s in full:
        r = {k: v for k, v in s["rep"]._asdict().items() if k != "Index"}
        r["appears_in"] = ", ".join(s["ranks"])
        rows.append(r)
    hl = pd.DataFrame(rows)
    if len(hl):
        show(f"Highlights: the same place listed once ({len(full)} places from "
             f"{sum(len(s['ranks']) for s in full)} leaderboard entries)", hl, lead="appears_in")

    # ---- files ------------------------------------------------------------------------
    def csv_ready(d, lead_cols):
        t = with_units(d)
        t["map"] = [maps_link(a, b) for a, b in zip(t.lat, t.lon)]
        t["mouth_map"] = [maps_link(a, b) for a, b in zip(t.mouth_lat, t.mouth_lon)]
        t = t.drop(columns=["cell_id", "mouth_id"])
        return t[lead_cols + [c for c in t.columns if c not in lead_cols]]

    parts = []
    for title, d in lists.items():
        if not len(d):
            continue
        t = d.copy()
        t.insert(0, "rank", np.arange(1, len(t) + 1))
        t.insert(0, "list", title)
        parts.append(csv_ready(t, ["list", "rank"]))
    if parts:
        pd.concat(parts).to_csv("leaderboards.csv", index=False)
        log("wrote leaderboards.csv")
    if len(hl):
        csv_ready(hl, ["appears_in"]).to_csv("highlights.csv", index=False)
        log("wrote highlights.csv")
    if KML_TOP > 0:
        write_kml("wrong_way_paths.kml", build_stories(lists, top=KML_TOP), fd, g)
    log("done")


if __name__ == "__main__":
    main()

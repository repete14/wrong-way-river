# Wrong-Way River

A game with real data: find the places on Earth where a raindrop would travel the furthest in the "wrong" direction. Think of a point 12 miles from the Pacific whose water runs 1,900 miles to the Bering Sea (the Yukon), or the Amazon's Andean headwaters, a short distance from the Pacific, that end up in the Atlantic.

For every land cell in a global flow-direction raster, the script computes:

- **flow**: how far the water actually travels to its endpoint.
- **straight**: the straight-line distance to the nearest ocean coast.
- **detour** = flow - straight (absolute "wrongness").
- **ratio** = flow / max(straight, 3 miles) (relative "wrongness").

It prints ranked lists, writes CSV files, and writes a KML file you can open in Google Earth to see the flow path and the straight line for the top results. Distances are in miles by default (kilometers with `WRONG_WAY_UNITS=km`).

## Quick start

1. **Install** (Python 3.9+; basic Python experience is enough):
   ```
   git clone <this repo>
   cd wrong-way-river
   pip install -e .
   ```
   This installs numpy, numba, scipy, rasterio, and pandas, and adds a `wrong-way-river` command. If you'd rather not install the project, `pip install -r requirements.txt` and run `python wrong_way_river.py` instead.
2. **Download the data** (next section) and put `hyd_glo_dir_30s.tif` in the folder you will run from.
3. **Run it:**
   ```
   wrong-way-river
   ```
   The first full global run took about 13 minutes on the author's machine. Reruns using the cache are much faster (see "Caching").
4. Read the console output, open `highlights.csv`, and open `wrong_way_paths.kml` in Google Earth.

`wrong-way-river --help` prints every setting, and `wrong-way-river --version` prints the version (currently 0.5.0).

Memory: peak RAM for the global file (about 725 million cells) hasn't been measured. The author's machine handled it, but if yours struggles, restrict the work with `WRONG_WAY_BBOX` (below). The cache folder uses roughly 6 GB for the global file, plus a scored-results file of a GB or two. It's safe to delete at any time.

## Data

You need the **HydroSHEDS v1 flow-direction (DIR) raster at 30 arc-seconds, global, as a GeoTIFF**.

- Start at the HydroSHEDS v1 product page: https://www.hydrosheds.org/products/hydrosheds. It links to the core data downloads page.
- Choose the flow-direction (DIR) layer, 30 arc-second resolution, global coverage, GeoTIFF format. The file inside the zip should be named `hyd_glo_dir_30s.tif`. Regional files (for example `hyd_na_dir_30s.tif` for North America) also work, but a basin that runs off the edge of a file cannot be scored. The North America file, for instance, cuts off the Yukon.
- Don't use the legacy archive (https://www.hydrosheds.org/downloads-archive, BIL/GRID formats) unless you have to, and don't use HydroSHEDS v2. v2 is 1 arc-second, about 900 times as many cells, and far too heavy for this script.
- HydroSHEDS data are free to use under the HydroSHEDS license; see the HydroSHEDS site for the terms and citation.

If the script says it can't find the flow-direction file, it hasn't found this download: put the file in the folder you run from, or set `WRONG_WAY_TIF` to its full path.

What the raster contains: eight direction codes (1, 2, 4, 8, 16, 32, 64, 128), 0 for outlets to the ocean, and 255 for nodata. The ocean is nodata, and, according to the Earth Engine catalog entry for this dataset, isolated 255 cells also mark inland sinks (endorheic basins).

## Settings

Everything is an environment variable, and everything is optional.

| Variable | Default | Meaning |
|---|---|---|
| `WRONG_WAY_TIF` | `hyd_glo_dir_30s.tif` | Path to the flow-direction GeoTIFF |
| `WRONG_WAY_BBOX` | whole file | `west,south,east,north` in degrees. Only this part of the file is read. It must contain whole basins **and** the coasts near them, because coasts outside the box are invisible and would make straight-line distances too large |
| `WRONG_WAY_REGION` | everywhere | `west,south,east,north`. The whole file is still used for flow paths and coasts, but only start cells inside this box are ranked. Use this for an "Asia only" or "Europe only" leaderboard. Doesn't handle boxes crossing 180 degrees longitude |
| `WRONG_WAY_UNITS` | `mi` | `mi` or `km`, for everything printed or written |
| `WRONG_WAY_TOP_N` | 20 | Results per leaderboard |
| `WRONG_WAY_PARETO_N` | 100 | Maximum entries in each Pareto frontier list |
| `WRONG_WAY_KML_TOP` | 10 | Entries per non-Pareto leaderboard drawn in the KML; `0` skips the KML. Pareto entries are always all drawn |
| `WRONG_WAY_FORCE_OCEAN` | none | `lat,lon;lat,lon`: treat the river mouth nearest each point as ocean-draining (see "Fixing misclassified basins") |
| `WRONG_WAY_CACHE_DIR` | `wrong_way_cache` | Cache folder |
| `WRONG_WAY_REFRESH` | unset | Set to anything to ignore the caches and recompute everything |

**Distance settings** can be given in miles or kilometers by adding `_MI` or `_KM` to the name (miles win if both are set). For example `WRONG_WAY_BASIN_SEP_MI=400` or `WRONG_WAY_BASIN_SEP_KM=650`.

| Name | Default | Meaning |
|---|---|---|
| `WRONG_WAY_DECLUSTER` | 50 km (31 mi) | Keep only the best result within this radius, so one ridge appears once |
| `WRONG_WAY_PARETO_DECLUSTER` | 15 km (9 mi) | Same, for the Pareto lists. Smaller, so frontier points near each other aren't lost |
| `WRONG_WAY_BASIN_SEP` | 1000 km (621 mi) | Two results at the same river mouth are both listed only if they are at least this far apart. `0` turns the basin rule off |
| `WRONG_WAY_MOUTH_MERGE` | 100 km (62 mi) | Big river mouths closer than this are treated as one, so a delta's many channels count as one basin |
| `WRONG_WAY_BAND` | 100 km (62 mi) | Near-coast band used by the speed-up (see "How it works") |
| `WRONG_WAY_COAST_CUTOFF` | 30 km (19 mi) | A flow path ending within this distance of the coast counts as ocean-draining; farther counts as endorheic |

Examples:

```powershell
# PowerShell: only rank starting points in Asia
$env:WRONG_WAY_REGION = "60,-12,150,55"
wrong-way-river
```

```
:: cmd
set WRONG_WAY_TIF=C:\data\hyd_glo_dir_30s.tif
set WRONG_WAY_BASIN_SEP_MI=0
wrong-way-river
```

```bash
# bash
WRONG_WAY_UNITS=km WRONG_WAY_BBOX="-170,14,-50,84" wrong-way-river
```

## What you get

**Console output**, in this order:
1. Progress with timings and cell counts.
2. A table of the biggest endorheic (closed-basin) sinks, as a sanity check. Real closed basins (Lake Chad, Lake Eyre, the Aral Sea and so on) should dominate it.
3. A "borderline" table: sinks classified endorheic but within 93 miles of the coast. These are the likeliest misclassifications.
4. Six leaderboards: ocean-draining and endorheic, each by absolute detour, by ratio, and as a Pareto frontier (results where no other result beats them on both measures). Each row has a `map` link to the location and a `mouth` link to where the water ends up.
5. An exactness report (see "How it works").
6. **Highlights**: the same places from all six lists merged, so each place is listed once with every list it appears in (for example `Pareto #3, Detour #1, Ratio #7`). Start here.

**`highlights.csv`**: the highlights table. **`leaderboards.csv`**: every row of every leaderboard, with river-mouth and nearest-coast coordinates.

**`wrong_way_paths.kml`**: one file with every highlight on it, for scanning in Google Earth.
- Highlights are grouped into three folders: the Pareto frontier (ocean-draining), other top ocean-draining results, and endorheic results. Each place is its own folder, so you can switch places on and off in the Places panel.
- Each place is labeled with its best rank, how many other lists it appeared in, and its numbers, for example `Pareto #3 (+2) · 1,551 mi · 6x`. The folder name and description list every rank it earned and link to each location.
- Colors show importance. The Pareto frontier is drawn bold, with a **red** flow path to the river mouth and a **green** straight line to the nearest coast. Other ocean-draining results are drawn thinner in muted orange and blue, and endorheic results in gray.
- All Pareto entries are drawn; the other lists draw their top 10 (`WRONG_WAY_KML_TOP`).
- Flow paths are thinned to about 800 vertices each to keep the file small.

### Opening the KML in Google Earth

The KML is a file, not a link. Either:
- **Google Earth Pro (desktop):** install it from https://www.google.com/earth/ if needed, then double-click `wrong_way_paths.kml`, or use File > Open.
- **Google Earth on the web:** use its Projects menu and look for the option to import a KML file. The exact wording may differ.

The map links in the console and CSV files open Google Maps directly and need no file.

## Caching

Everything slow is cached per input file and bounding box in `wrong_way_cache/`:

1. **Flow results** (reading the raster and walking every flow path). This is the slowest step and is done once.
2. **River mouth table.**
3. **Scored cells**, which depend on the cutoff, forced-ocean list, band, region, and minimum thresholds.

Changing only the grouping and display settings (`TOP_N`, `DECLUSTER`, `BASIN_SEP`, `MOUTH_MERGE`, `KML_TOP`, `UNITS`) normally reuses the scored cells and finishes in seconds. If your new grouping settings are stricter than the cached run assumed, so cells the cache skipped could now matter, the script says so and rescores automatically.

## Fixing misclassified basins

A flow path that ends within 19 miles (30 km) of the coast counts as ocean-draining. Big deltas and lakes with outlets can end farther inland than that and then show up as "endorheic". Use the borderline table to spot them.

- Raise `WRONG_WAY_COAST_CUTOFF`, which affects every basin. A very large cutoff mislabels real closed basins that sit near coasts (for example, the Dead Sea).
- Or force specific mouths: `WRONG_WAY_FORCE_OCEAN="lat,lon;lat,lon"`, using the lat/lon printed in the borderline table. The mouth nearest each point (within 15 miles) is treated as ocean-draining.

The tables also show a `HydroSHEDS_flag` column: `outlet(0)` means the dataset itself marks that cell as an ocean outlet, and `inland-sink` means it marks a sink. That is useful evidence when deciding whether a borderline basin is genuine.

Example from a global run: the largest basins flagged endorheic included what looked like the Murray-Darling terminal (near -35.35, 139.38) and the Irrawaddy delta (near 16.36, 95.04), both of which probably reach the sea. Check these in the borderline table before treating them as closed basins.

## How it works

This section explains the ideas in enough detail to rebuild the script from scratch. Each step maps to a function in `wrong_way_river.py`.

### 1. Turn the question into two numbers per cell

HydroSHEDS tells us, for every ~1 km land cell, which of its 8 neighbors water flows to (a "D8" flow-direction grid). Treat every cell as a raindrop and measure:

- the **flow distance**: how far the drop travels following those arrows until it reaches an endpoint, and
- the **straight-line distance** from the starting cell to the nearest ocean coast.

The "wrong way" gap is the difference (detour) or the quotient (ratio). Doing this per cell means the script never has to decide where a river "starts", which is a fuzzy question.

### 2. Flow distance: each cell is solved once (`flow_distance`)

The grid is really a forest of trees: every cell has one downstream neighbor, and the roots are river mouths and inland sinks. A cell's flow distance is the step to its downstream neighbor plus that neighbor's flow distance, so the cells nearest the roots must be solved first.

Instead of ordering the whole grid, the code does memoized path-following. For each unsolved cell it walks downstream, pushing cells on a stack, until it reaches a cell that is already solved, a sink, the ocean, or the edge of the region. Then it unwinds the stack back toward the start, filling in each cell's distance and recording the endpoint cell ("terminal cell") its water ends in. Every cell is written exactly once, so the walk is linear in the number of cells. It's compiled with numba because plain Python would take hours.

Details: the length of a step depends on latitude (cells get narrower toward the poles, so east-west steps are scaled by the cosine of latitude, and diagonals use the hypotenuse). Paths that leave the bounding box are marked "truncated" and excluded, because their true endpoint is unknown.

### 3. Straight-line distance: a KD-tree on a sphere (`straight_for`)

- **Finding the ocean:** every connected blob of nodata cells larger than 5,000 cells is ocean (connected-component labeling). Small nodata specks are mostly the dataset's inland-sink markers, not ocean.
- **Finding the coast:** coast cells are land cells touching the ocean (one step of binary dilation).
- **Measuring:** every coast cell becomes a 3D point on a sphere of radius 6,371 km, and a KD-tree answers "nearest coast point" queries. Doing this in 3D instead of on the latitude/longitude grid avoids the distortion near the poles. The tree returns a straight-through-the-earth chord, which is converted to a true great-circle distance.

### 4. Ranking: absolute, relative, and the Pareto frontier

Absolute detour favors huge systems, such as a 3,000-mile river that starts 25 miles from the sea. Ratio favors tiny absurdities, such as a stream 2 miles from the coast that runs 200 miles. Neither is "right", so the script reports both, plus the **Pareto frontier**.

A cell is on the frontier if no other cell beats it on both measures at once. Picture a scatter plot with detour on one axis and ratio on the other: the frontier is the outer edge of the cloud. Anything on it is the "most wrong" cell under some weighting of the two measures, and the ends of the frontier are the interesting characters. To find it, sort cells by detour, largest first, and keep each cell whose ratio is higher than that of every cell before it.

The ratio uses a floor of 5 km (3 miles) on the straight-line distance, so cells right on the coast (distance near zero) don't produce infinite ratios.

### 5. Ocean-draining versus endorheic

Water that ends in an inland basin never reaches an ocean, but a raindrop near a coast that ends in a desert lake is still very wrong. So both kinds are scored the same way, against the straight line to the ocean, and listed separately. The classification is made per river mouth or sink: if the terminal cell is within 30 km of the coast, the whole basin counts as ocean-draining, otherwise as endorheic. This is a heuristic. Big deltas and lakes with outlets can fall on the wrong side, which is why the script prints a "borderline" table and lets you force mouths with `WRONG_WAY_FORCE_OCEAN`. The dataset's own flag (0 for ocean outlet, inland-sink marker otherwise) is shown beside each mouth as a second opinion.

### 6. The speed-up: most cells can be skipped (`score_cells`)

A KD-tree query for every one of ~170 million candidate cells is the slow part. The key fact is that detour = flow - straight and straight is never negative, so **a cell's detour can never exceed its flow distance**, and the flow distances are already known.

- **Phase 1:** query only cells near the coast. The tree is given a distance limit (`WRONG_WAY_BAND`, default 62 miles) and bails out early for cells beyond it. This finds all the high-ratio candidates, since a high ratio needs a short straight-line distance.
- **Phase 2:** take the N-th best detour found so far, call it T. Any inland cell whose flow distance is below T cannot possibly make the list, so only inland cells with flow above T get a full query.

In the first global run, about 170 million cells had flows over 62 miles and only about 34 million (20%) needed a straight-line distance. The absolute-detour lists are exact under this scheme. It relies on the N-th best value only ever going up as more cells are scored, so a threshold from a partial set is always safe. That is true for the plain spatial grouping; I believe it also holds with the basin rule, but I haven't proved it.

For ratio lists, any skipped cell is farther than the band from the coast, so its ratio is below flow / band. If the N-th best ratio exceeds that bound, the list is provably exact. The report at the end checks this and tells you what band to use if it isn't.

A side fact that the code doesn't need: following water upstream never decreases the absolute detour. If cell U drains into cell X, U's flow distance grows by the step, but its straight-line distance can grow by at most that step (the triangle inequality). So the very top absolute results are always tips of tributaries. A ratio can't be treated the same way, because moving upstream often means moving away from the coast.

### 7. Grouping results (`decluster`, `build_mouth_groups`, `build_stories`)

Three layers of grouping keep the lists from showing the same thing over and over.

- **Distance:** the best cell's neighbors are almost as wrong, so a raw top 20 would be one hillside twenty times. The script takes the best cell, discards everything within 31 miles, and repeats.
- **Basin rule:** each cell records its terminal cell, so we know which river mouth it belongs to. Two results with the same mouth are both kept only if they are at least 621 miles apart. That merges near-duplicates (the Nile basin has several spots a few hundred miles apart) but still lets two very different routes into one basin appear. Deltas have many channel mouths, so big mouths (those draining at least 20,000 cells) closer than 62 miles are first merged into one basin group by connected components. Without that step, the Nile or Rio de la Plata would count as several basins.
- **Highlights:** the same place can rank in several lists (for example first by detour and third on the Pareto frontier). Entries are merged into one "story" when they are within the distance radius of each other, or share a basin group and are within the basin separation. Lists are visited in order of importance (ocean Pareto, detour, ratio, then the endorheic lists), and the first entry of a story supplies the cell that is drawn.

The Pareto lists use a smaller merge radius (9 miles) and no basin rule, since the goal there is the full shape of the frontier. They are still merged into stories in the KML and highlights.

### 8. Drawing the paths (`trace_path`, `write_kml`)

The direction raster isn't kept in the cache, only the flow distances. That is enough: standing on a cell, the downstream neighbor is the one whose flow distance plus the step to reach it equals the current cell's distance. Repeat until the distance reaches zero (the mouth). The path is thinned to about 800 points, and the straight line to the nearest coast point is drawn alongside it. If a trace fails to reach the mouth, the KML description says so and the script logs a count.

### 9. Caching

Flow results, the table of river mouths, and the scored cells are saved to `wrong_way_cache/`. The scored cells are tagged with the settings that affect them (cutoff, forced mouths, band, region). When you change only grouping settings, the script checks that the cache still contains every cell those settings could need, and reuses it if so.

## Caveats

- Resolution is 30 arc-seconds (about 1 km), which is fine for the game but not for exact headwater locations.
- HydroSHEDS v1 is built from SRTM elevation data, which stops at 60 N. Above that, the data come from a coarser, lower-quality source, so high-latitude results are less reliable.
- The coast comes from the raster's own nodata mask, so fjords and inlets count as coast.
- Paths that leave the bounding box can't be scored and are dropped, so use a box that holds whole basins and their nearby coasts.
- Flow paths follow the DEM. Real-world quirks like canals, dams, and groundwater are not represented.
- The code was tested end to end on small synthetic rasters and, in earlier versions, on the real global file. The mouth merging, highlights, color scheme, region filter, and miles/kilometers settings have only been exercised on synthetic data.

## License

Licensed under the [PolyForm Noncommercial License 1.0.0](https://polyformproject.org/licenses/noncommercial/1.0.0) (SPDX: `PolyForm-Noncommercial-1.0.0`). In plain English:

- You can use it, copy it, share it, fork it, and modify it for any **noncommercial** purpose: personal use, hobby projects, research, and education.
- You can't use it, or a modified version of it, for commercial purposes (for example, selling it or building it into something you sell) without separate permission from the author.
- If you pass it on, include the license (or a link to it) and the copyright notice in `LICENSE`.

This summary is not the license; the linked license text is the authority, and it is not legal advice. The HydroSHEDS data you download is covered by its own license (see "Data" above), which is separate from this one.

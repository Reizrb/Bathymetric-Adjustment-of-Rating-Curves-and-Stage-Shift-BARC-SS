# BARC-SS: Bathymetric Adjustment of Rating Curves and Stage Shift

Python code for **BARC** and **BARC-SS**, which add the river channel geometry missing
from DEMs to the synthetic rating curves (SRCs) of the
[OWP HAND-FIM](https://github.com/NOAA-OWP/inundation-mapping) flood inundation
mapping framework.

DEMs represent the water surface, not the channel bed, so HAND-derived rating curves
miss the in-channel cross-sectional area and wetted perimeter and underestimate
channel capacity. NOAA OWP's
[Bathymetric Adjustment of Rating Curves (BARC)](https://github.com/NOAA-OWP/inundation-mapping/blob/dev/src/bathymetric_adjustment.py)
module corrects the rating curves by adding channel geometry to the hydrotable.
**BARC-SS** keeps that channel-geometry adjustment and adds a reach-specific stage
shift based on the estimated mean-flow depth.

Both methods here use machine-learning estimates of bankfull width, bankfull depth,
and mean-flow depth for every NHDPlusV2.1 reach across CONUS
([HydroGeom](https://github.com/Reizrb/HydroGeom); Zarrabi et al., 2025), with a
channel shape you choose.

| Method | What it does |
|--------|--------------|
| **BARC** | Adds the missing bankfull cross-sectional area and wetted perimeter to each reach's hydrotable and recomputes discharge with Manning's equation |
| **BARC-SS** | BARC, plus shifts the stage of each reach's rating curve up by its mean-flow depth and re-anchors the curve at zero flow |

## How it works

BARC-SS has five steps. BARC is Steps 1 to 3.

### Step 1. Generate the initial hydrotable

The initial hydrotable of each reach comes from OWP HAND-FIM version 4.5, through
[FIMserv](https://github.com/sdmlua/FIMserv). It holds the geometric and hydraulic
attributes used to build the rating curve at each stage (S): cross-sectional flow area
(A), volume (V), top width (B), bed area (BA), wetted perimeter (P), hydraulic radius
(R), and discharge (Q). Because the DEM captures only the channel above the water
surface (depth d<sub>1</sub>), the initial rating curve underestimates the true
channel capacity.

<p align="center"><img src="assets/fig4a_step1_initial.png" alt="Step 1: the initial HAND-derived rating curve captures only the DEM-visible channel of depth d1" width="90%"></p>
<p align="center"><sub><b>Figure 1a.</b> The initial HAND-derived SRC (red) captures only the DEM-visible geometry of depth d<sub>1</sub>, underestimating the true channel capacity.</sub></p>

### Step 2. Apply the BARC channel-geometry adjustment

For the chosen channel shape, the bankfull cross-sectional area A<sub>bnk</sub> and
wetted perimeter P<sub>bnk</sub> are computed from each reach's bankfull width
w<sub>bnk</sub> and depth d<sub>bnk</sub> with HydroGeomKit's `xsec_area` and
`wetted_perimeter` (see [Channel shapes](#channel-shapes)). This in-channel geometry
is missing from the DEM (A<sub>missing</sub> = A<sub>bnk</sub>,
P<sub>missing</sub> = P<sub>bnk</sub>) and is added to the hydrotable at every stage:

| | |
|---|---|
| V<sub>adjusted</sub> = V<sub>initial</sub> + A<sub>missing</sub> × L × 1000 | (1) |
| BA<sub>adjusted</sub> = BA<sub>initial</sub> + P<sub>missing</sub> × L × 1000 | (2) |
| P<sub>adjusted</sub> = P<sub>initial</sub> + P<sub>missing</sub> | (3) |
| A<sub>adjusted</sub> = A<sub>initial</sub> + A<sub>missing</sub> | (4) |

where L is the reach length (km) and 1000 converts it to meters.

### Step 3. Recompute the BARC rating curve

Hydraulic radius and discharge are recomputed from the adjusted area and wetted
perimeter:

| | |
|---|---|
| R<sub>adjusted</sub> = A<sub>adjusted</sub> / P<sub>adjusted</sub> | (5) |
| Q<sub>adjusted</sub> = (1/n) × A<sub>adjusted</sub> × R<sub>adjusted</sub><sup>2/3</sup> × S<sub>0</sub><sup>1/2</sup> | (6) |

where n is Manning's roughness coefficient and S<sub>0</sub> the channel slope. Both are
kept from the initial hydrotable, so the change in the rating curve reflects only the
channel geometry. The BARC rating curve is anchored at stage = 0 and discharge = 0.
**This completes BARC.**

<p align="center"><img src="assets/fig4b_step2-3_barc.png" alt="Steps 2 and 3: the missing area and perimeter are added and discharge is recomputed, shifting the curve right by Q missing" width="90%"></p>
<p align="center"><sub><b>Figure 1b.</b> Bankfull geometry gives A<sub>missing</sub> and P<sub>missing</sub>, which are added to the hydrotable (Eqs. 1–4). R and Q are then recomputed with Manning's equation (Eqs. 5–6), shifting the BARC curve (yellow) right by Q<sub>missing</sub>.</sub></p>

### Step 4. Apply the mean-flow-depth stage adjustment (BARC-SS)

The estimated mean-flow depth d<sub>mf</sub> is added to the stage values of each
reach's BARC rating curve:

| | |
|---|---|
| Stage<sub>BARC-SS</sub> = Stage<sub>BARC</sub> + d<sub>mf</sub> | (7) |

Using d<sub>mf</sub> does not assume the elevation data were collected under mean-flow
conditions. It provides a spatially continuous, reach-specific estimate that sets the
size of the stage adjustment. The actual river stage when the elevation data were
collected may differ, and because d<sub>mf</sub> is estimated for each reach, the
adjustment varies across the river network.

### Step 5. Anchor the BARC-SS rating curve at zero flow

The stage shift moves the whole BARC rating curve up, including its zero-flow point.
A new point with stage = 0 and discharge = 0 is therefore added at the start of each
BARC-SS rating curve, so it starts at the origin and is continuous from zero flow.
Because this point is joined to the first shifted stage, BARC-SS rating curves have a
nearly straight segment at very low discharges.

<p align="center"><img src="assets/fig4c_step4-5_barc-ss.png" alt="Steps 4 and 5: stage is shifted by the mean-flow depth and the curve is anchored at zero flow" width="90%"></p>
<p align="center"><sub><b>Figure 1c.</b> Stage values are shifted by the mean-flow depth d<sub>mf</sub> to reference the channel bed, and a zero-flow anchor (S = 0, Q = 0) is added, giving the final BARC-SS curve (yellow).</sub></p>

### Flood maps

The adjusted hydrotables are used in OWP HAND-FIM as usual, and the HAND raster is
never modified. For each reach, the stage for the forcing discharge is linearly
interpolated from the rating curve and compared with HAND to find the inundated area.

### Channel shapes

The channel shape determines how bankfull width and depth become cross-sectional area
and wetted perimeter. All cross sections are symmetric about the channel centerline.

- **Trapezoid, parameter a** (`shape="a"`): a = w<sub>bot</sub> / w<sub>bnk</sub>, the
  bottom width divided by the bankfull width. a = 0 is triangular, a = 1 rectangular,
  and values in between trapezoidal. A<sub>bnk</sub> and P<sub>bnk</sub> follow the
  standard trapezoid relationships with w<sub>bot</sub> = a × w<sub>bnk</sub>.
- **Dingman's r** (`shape="r"`; Dingman & Afshari, 2018): bank curvature follows
  z = d<sub>bnk</sub> (2/w<sub>bnk</sub>)<sup>r</sup> x<sup>r</sup> for
  0 ≤ x ≤ w<sub>bnk</sub>/2. r = 1 is triangular, r = 2 parabolic, and larger r
  approaches rectangular. A<sub>bnk</sub> = [r / (r + 1)] × w<sub>bnk</sub> × d<sub>bnk</sub>,
  and P<sub>bnk</sub> is the arc length of the profile, computed by adaptive numerical
  integration.
- **Side slope z** (`shape="z"`): banks of slope z:1, also available in HydroGeomKit.

The two parameters give the same cross-sectional area when a = (r − 1) / (r + 1).

The paper tests these scenarios:

| Parameter a | S01 | S02 | S03 | S04 | S05 | S06 | |
|---|---|---|---|---|---|---|---|
| Value | 0.0 | 0.2 | 0.4 | 0.6 | 0.8 | 1.0 | |
| **Parameter r** | **R01** | **R02** | **R03** | **R04** | **R05** | **R06** | **R07** |
| Value | 1 | 2 | 3 | 4 | 6 | 10 | 100 |

No reach-specific channel-shape data exist across CONUS, so these scenarios span
idealized shapes from triangular to rectangular to measure how much the choice of shape
affects flood maps, not to pick a single best shape.

### Results in brief

Across six HUC8 case studies, the paper finds:

- Larger channel cross-sectional area reduces predicted flood extent. Flood maps are most
  sensitive to channel shape for frequent floods (2- to 10-year), and often for lower- to
  intermediate-order streams.
- BARC alone generally reduces flood extent and can increase underprediction.
- BARC-SS improves agreement with remote-sensing benchmark maps for the Arkansas,
  Illinois, and Neuse Rivers, raising mean CSI from 0.32, 0.30, and 0.39 to 0.63, 0.64,
  and 0.74. In the other three case studies, better flood detection comes with more false
  alarms.
- BARC-SS is less sensitive to the channel-shape assumption than BARC: the mean CSI range
  across the six a scenarios drops from 0.20 to 0.047.
- No single channel shape or adjustment performs best everywhere.

## Install

```bash
pip install "git+https://github.com/Reizrb/Bathymetric-Adjustment-of-Rating-Curves-and-Stage-Shift-BARC-SS"
```

To run the example notebook, which also needs [FIMserv](https://github.com/sdmlua/FIMserv)
and Jupyter:

```bash
pip install "barcss[example] @ git+https://github.com/Reizrb/Bathymetric-Adjustment-of-Rating-Curves-and-Stage-Shift-BARC-SS"
```

Python 3.10 to 3.13.

## Quick start

```python
import fimserve as fm
import barcss

huc = "11110203"
fm.DownloadHUC8(huc, version="4.5")                     # OWP HAND-FIM outputs into ./output

ids = barcss.feature_ids("output", huc)                  # reaches in the rating curves
geometry = barcss.bankfull_geometry(ids, shape="a", a=0.4)
barcss.apply("output", huc, geometry, method="BARC-SS")  # or method="BARC"

fm.getNWMretrospectivedata(huc, "2015-12-26", "2016-01-05", value_times=["2016-01-03 00:00:00"])
fm.runOWPHANDFIM(huc)
barcss.plot.flood_map("output", huc)                     # map of the inundated area
```

## Example

[`examples/BARC-SS_HUC8_11110203_a0.4.ipynb`](examples/BARC-SS_HUC8_11110203_a0.4.ipynb)
reproduces one configuration from the paper: **Case Study 1, the Arkansas River
(HUC8 11110203)**, during the 2015–2016 Midwest U.S. flooding, forced by hourly NWM
retrospective streamflow on **3 January 2016 at 00:00**, with a trapezoidal channel of
**a = 0.4 (scenario S03)**. It downloads HAND-FIM, gets the channel geometry, applies
BARC-SS, plots a rating curve before and after, runs HAND-FIM, and maps the flood. Run
it from an empty folder, since FIMserv writes `code/`, `data/`, and `output/` folders
where it runs.

A benchmark flood map for this event is available in
[FIMbench](https://fimbench.tethys.ciroh.org/) to evaluate the result.

## Functions

| Function | What it does |
|----------|--------------|
| `feature_ids(output_dir, huc8)` | NWM feature IDs (COMIDs) in the HUC8's rating curves |
| `bankfull_geometry(comids, shape, **param)` | Bankfull and mean-flow geometry from HydroGeom, plus A<sub>bnk</sub> and P<sub>bnk</sub> for the shape |
| `apply(output_dir, huc8, geometry, method)` | Apply BARC or BARC-SS to every branch of the HUC8 |
| `restore(output_dir, huc8)` | Put back the original HAND-FIM files |
| `keep_flood_map(output_dir, huc8, label)` | Rename the latest flood map so FIM can run again with other settings |
| `plot.rating_curve(output_dir, huc8, hydroid)` | Original and adjusted rating curve of one reach |
| `plot.flood_map(output_dir, huc8)` | Map of the inundated cells with the HUC8 boundary |

`apply` saves the untouched HAND-FIM files in `output/flood_<huc8>/barcss_original/`
the first time and always starts from them, so running it again (with the same or
different settings) never adds the geometry twice.

### Notes

- Reaches without a geometry prediction are left unadjusted (no added area, no stage shift).
- The geometry is added at every stage, including stage 0, so the BARC rating curve has
  two points at stage 0: the zero-flow anchor and the adjusted stage-0 discharge
  (the horizontal segment of length Q<sub>missing</sub> in Figure 1b).
- As in the OWP bathymetry adjustment, the subdivision columns (`subdiv_applied`,
  `overbank_n`, `channel_n`, `subdiv_discharge_cms`) are removed from the hydroTable,
  and the original discharge is kept as `discharge_cms_org`.

## Tests

```bash
pip install -e ".[test]"
python -m pytest tests
```

The tests check Equations 1 to 7, that BARC and BARC-SS differ only by the stage
shift, that re-running and restoring work, and that the results match the original
code used in the paper (`tests/reference/original_utils.py`) exactly.

## Citation

If you use this code, please cite:

- Zarrabi, R., Cohen, S., Pruitt, C., Baruah, A., McDermott, R., & Chen, Y. (2026).
  Sensitivity of terrain-based flood inundation model (OWP HAND-FIM) predictions to
  channel geometry: Insights from bathymetric adjustments of rating curves and stage
  shift. *Journal of Hydrology* (under review).

and, for the channel geometry:

- Zarrabi, R., McDermott, R., Erfani, S. M. H., & Cohen, S. (2025). Bankfull and
  mean-flow channel geometry estimation through machine learning algorithms across
  the CONtiguous United States (CONUS). *Water Resources Research*, 61(2).
  https://doi.org/10.1029/2024WR037997

Other references:

- Dingman, S. L., & Afshari, S. (2018). Field verification of analytical at-a-station
  hydraulic-geometry relations. *Journal of Hydrology*.

## Related repositories

| Repository | What it does |
|------------|--------------|
| [HydroGeom](https://github.com/Reizrb/HydroGeom) | Channel geometry for every NHDPlusV2.1 reach: web page, API, and Python package |
| [Bankfull-and-mean-flow-channel-geometry-for-CONUS](https://github.com/Reizrb/Bankfull-and-mean-flow-channel-geometry-for-CONUS) | The ML models and dataset behind the channel geometry |
| [FIMserv](https://github.com/sdmlua/FIMserv) | Runs OWP HAND-FIM for any HUC8 |

## License

[MIT](LICENSE)

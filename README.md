# BARC-SS: Bathymetric Adjustment of Rating Curves and Stage Shift

Python code to apply **BARC** and **BARC-SS** to the synthetic rating curves of
[OWP HAND-FIM](https://github.com/NOAA-OWP/inundation-mapping) and map the resulting floods.

- **BARC** adds the bankfull channel area and wetted perimeter missing from the DEM to each
  reach's hydrotable and recomputes discharge with Manning's equation.
- **BARC-SS** does the same, then shifts each rating curve's stage up by the reach's
  mean-flow depth and re-anchors it at zero flow.

Channel width and depth come from [HydroGeom](https://github.com/Reizrb/HydroGeom), and
you choose the channel shape used to turn them into area and wetted perimeter.

<p align="center"><img src="assets/fig4a_step1_initial.png" alt="Initial HAND-FIM rating curve" width="75%"><br><sub>Initial HAND-FIM: the DEM captures only the channel above the water surface.</sub></p>
<p align="center"><img src="assets/fig4b_step2-3_barc.png" alt="BARC: missing geometry added, discharge recomputed" width="75%"><br><sub>BARC: the missing channel geometry is added and discharge recomputed.</sub></p>
<p align="center"><img src="assets/fig4c_step4-5_barc-ss.png" alt="BARC-SS: stage shifted by mean-flow depth and anchored at zero flow" width="75%"><br><sub>BARC-SS: stage is also shifted by the mean-flow depth and anchored at zero flow.</sub></p>

## Requirements

- Python 3.10 to 3.13
- Internet access: HAND-FIM data are downloaded from the public `ciroh-owp-hand-fim` AWS
  bucket (no account needed), NWM discharge from the NWM retrospective, and channel
  geometry from the HydroGeom API
- A few GB of free disk space per HUC8

## Install

```bash
pip install "barcss[example] @ git+https://github.com/Reizrb/Bathymetric-Adjustment-of-Rating-Curves-and-Stage-Shift-BARC-SS"
```

This installs `barcss`, [HydroGeomKit](https://github.com/Reizrb/HydroGeom),
[FIMserv](https://github.com/sdmlua/FIMserv) (which runs HAND-FIM), and Jupyter.
Without `[example]`, only `barcss` and HydroGeomKit are installed.

## Reproduce the example

The example applies BARC-SS to **HUC8 11110203** (Arkansas River) with a trapezoidal
channel **a = 0.4**, and maps the flood of **3 January 2016, 00:00** using NWM
retrospective discharge and OWP HAND-FIM version 4.5.

```bash
mkdir barcss-example && cd barcss-example
curl -LO https://raw.githubusercontent.com/Reizrb/Bathymetric-Adjustment-of-Rating-Curves-and-Stage-Shift-BARC-SS/main/examples/BARC-SS_HUC8_11110203_a0.4.ipynb
jupyter notebook BARC-SS_HUC8_11110203_a0.4.ipynb
```

Run it from an empty folder: FIMserv creates `code/`, `data/`, and `output/` folders
where it runs. The notebook produces:

| Output | Where |
|--------|-------|
| Adjusted hydrotables | `output/flood_11110203/11110203/branches/*/hydroTable_*.csv` and `output/flood_11110203/11110203/hydrotable.csv` |
| Original (unadjusted) files | `output/flood_11110203/barcss_original/` |
| Flood map (GeoTIFF, EPSG:5070, 1 = flooded) | `output/flood_11110203/11110203_inundation/` |
| Flood map figure | `flood_map_11110203_BARC-SS_a0.4.png` |

## Use it for any HUC8, method, or channel shape

```python
import fimserve as fm
import barcss

huc = "11110203"
fm.DownloadHUC8(huc, version="4.5")                        # OWP HAND-FIM outputs into ./output

ids = barcss.feature_ids("output", huc)                     # reaches in the rating curves
geometry = barcss.bankfull_geometry(ids, shape="a", a=0.4)  # channel geometry for the chosen shape
barcss.apply("output", huc, geometry, method="BARC-SS")     # or method="BARC"

fm.getNWMretrospectivedata(huc, "2015-12-26", "2016-01-05", value_times=["2016-01-03 00:00:00"])
fm.runOWPHANDFIM(huc)
barcss.plot.flood_map("output", huc)
```

**Method:** `"BARC"` or `"BARC-SS"`.

**Channel shape** (passed to HydroGeomKit):

| Shape | Parameter | Range |
|-------|-----------|-------|
| `shape="a"` | trapezoid, bottom width / top width | `a=0` triangle to `a=1` rectangle |
| `shape="r"` | Dingman's shape parameter | `r=1` triangle, `r=2` parabola, larger `r` toward rectangle |
| `shape="z"` | side slope z:1 | `z=0` rectangle |

**Running again with other settings.** `apply` saves the original HAND-FIM files the first
time and always starts from them, so you can rerun with another method or shape without
downloading again. FIMserv skips a run if a flood map for the same discharge already
exists, so keep the previous map under its own name first:

```python
barcss.keep_flood_map("output", huc, "BARC-SS_a0.4")
geometry = barcss.bankfull_geometry(ids, shape="a", a=1.0)
barcss.apply("output", huc, geometry, method="BARC")
fm.runOWPHANDFIM(huc)
```

To map the original HAND-FIM for comparison, run `barcss.restore("output", huc)` and then
`fm.runOWPHANDFIM(huc)`.

## Functions

| Function | What it does |
|----------|--------------|
| `feature_ids(output_dir, huc8)` | NWM feature IDs (NHDPlusV2.1 COMIDs) in the HUC8's rating curves |
| `bankfull_geometry(comids, shape, **param)` | Bankfull and mean-flow width and depth from HydroGeom, plus bankfull area and wetted perimeter for the shape |
| `apply(output_dir, huc8, geometry, method)` | Apply BARC or BARC-SS to every branch of the HUC8 |
| `restore(output_dir, huc8)` | Put back the original HAND-FIM files |
| `keep_flood_map(output_dir, huc8, label)` | Rename the latest flood map so FIM can run again |
| `plot.rating_curve(output_dir, huc8, hydroid)` | Original and adjusted rating curve of one reach |
| `plot.flood_map(output_dir, huc8)` | Map of the flooded cells with the HUC8 boundary |

## What the code does

For each reach, with bankfull area A<sub>bnk</sub> and wetted perimeter P<sub>bnk</sub>
from the chosen shape, reach length L (km), and the hydrotable's Manning's n and slope S
(both unchanged):

| | |
|---|---|
| Volume += A<sub>bnk</sub> × L × 1000 | (1) |
| Bed area += P<sub>bnk</sub> × L × 1000 | (2) |
| Wetted perimeter P += P<sub>bnk</sub> | (3) |
| Wet area A += A<sub>bnk</sub> | (4) |
| R = A / P | (5) |
| Q = (1/n) × A × R<sup>2/3</sup> × S<sup>1/2</sup> | (6) |
| BARC-SS only: stage += mean-flow depth d<sub>mf</sub> | (7) |

Both methods then add a zero-flow point (stage = 0, discharge = 0) at the start of each
rating curve. The HAND raster is not changed.

- Reaches without a geometry prediction are left unadjusted.
- The subdivision columns (`subdiv_applied`, `overbank_n`, `channel_n`,
  `subdiv_discharge_cms`) are removed from the hydroTable, as in the OWP bathymetry
  adjustment, and the original discharge is kept as `discharge_cms_org`.

## Tests

```bash
git clone https://github.com/Reizrb/Bathymetric-Adjustment-of-Rating-Curves-and-Stage-Shift-BARC-SS
cd Bathymetric-Adjustment-of-Rating-Curves-and-Stage-Shift-BARC-SS
pip install -e ".[test]"
python -m pytest tests
```

The tests run on a small built-in HUC8 (no downloads). They check Equations 1 to 7, that
BARC and BARC-SS differ only by the stage shift, that reruns and `restore` work, and that
the results match the original research code (`tests/reference/original_utils.py`).

## Citation

If you use this code, please cite:

- Zarrabi, R., Cohen, S., Pruitt, C., Baruah, A., McDermott, R., & Chen, Y. (2026).
  Sensitivity of terrain-based flood inundation model (OWP HAND-FIM) predictions to
  channel geometry: Insights from bathymetric adjustments of rating curves and stage
  shift. *Journal of Hydrology* (under review).
- Zarrabi, R., McDermott, R., Erfani, S. M. H., & Cohen, S. (2025). Bankfull and
  mean-flow channel geometry estimation through machine learning algorithms across
  the CONtiguous United States (CONUS). *Water Resources Research*, 61(2).
  https://doi.org/10.1029/2024WR037997

## Related repositories

| Repository | What it does |
|------------|--------------|
| [HydroGeom](https://github.com/Reizrb/HydroGeom) | Channel geometry for every NHDPlusV2.1 reach: web page, API, and Python package |
| [FIMserv](https://github.com/sdmlua/FIMserv) | Downloads and runs OWP HAND-FIM for any HUC8 |
| [OWP inundation-mapping](https://github.com/NOAA-OWP/inundation-mapping) | OWP HAND-FIM, including its bathymetric adjustment module |

## License

[MIT](LICENSE)

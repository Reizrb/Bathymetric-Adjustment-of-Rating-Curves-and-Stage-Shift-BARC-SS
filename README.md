# BARC-SS: Bathymetric Adjustment of Rating Curves and Stage Shift

Python code for **BARC** and **BARC-SS**, two methods that add the river channel
geometry missing from DEMs to the synthetic rating curves (SRCs) of the
[OWP HAND-FIM](https://github.com/NOAA-OWP/inundation-mapping) flood inundation
mapping framework.

DEMs capture the water surface, not the channel bed, so HAND-FIM rating curves miss
the in-channel cross-sectional area and wetted perimeter. This underestimates channel
conveyance and can overestimate flood extent. BARC and BARC-SS add that missing
geometry, using machine-learning predictions of bankfull and mean-flow channel width
and depth for every NHDPlusV2.1 reach across CONUS
([HydroGeom](https://github.com/Reizrb/HydroGeom); Zarrabi et al., 2025).

| Method | What it does |
|--------|--------------|
| **BARC** | Adds the missing bankfull area and wetted perimeter to each reach's hydrotable and recomputes discharge with Manning's equation |
| **BARC-SS** | BARC, plus shifts the stage datum up by the mean-flow depth so stage is measured from the channel bed rather than the DEM water surface |

## How it works

For each reach and channel shape, the bankfull cross-sectional area A<sub>bnk</sub>
and wetted perimeter P<sub>bnk</sub> are computed from the predicted bankfull width
and depth. These are the geometry missing from the DEM
(A<sub>missing</sub> = A<sub>bnk</sub>, P<sub>missing</sub> = P<sub>bnk</sub>) and are
added to every stage of the initial hydrotable:

| | |
|---|---|
| V<sub>adj</sub> = V + A<sub>missing</sub> · L · 1000 | (1) |
| B<sub>adj</sub> = B + P<sub>missing</sub> · L · 1000 | (2) |
| P<sub>adj</sub> = P + P<sub>missing</sub> | (3) |
| A<sub>adj</sub> = A + A<sub>missing</sub> | (4) |
| R<sub>adj</sub> = A<sub>adj</sub> / P<sub>adj</sub> | (5) |
| Q<sub>adj</sub> = (1/n) · A<sub>adj</sub> · R<sub>adj</sub><sup>2/3</sup> · S<sup>1/2</sup> | (6) |

where V is volume, B bed area, L reach length (km), n Manning's roughness, and S
slope (n and S are unchanged). **BARC-SS** then shifts every stage up by the
reach's mean-flow depth d<sub>mf</sub>. Both methods add a zero-flow anchor
(stage = 0, discharge = 0) at the start of each rating curve.

### Channel shapes

The channel shape sets how bankfull width and depth become area and wetted
perimeter. Any shape supported by
[HydroGeomKit](https://github.com/Reizrb/HydroGeom#channel-shape-parameters) can be used:

- `shape="a"`: trapezoid with bottom width = a × top width (`a=0` triangle, `a=1` rectangle)
- `shape="r"`: Dingman's shape parameter (`r=1` triangle, `r=2` parabola)
- `shape="z"`: side slope z:1

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
runs the full workflow for HUC8 **11110203** (Arkansas River basin) with a trapezoidal
channel (**a = 0.4**) and the 3 January 2016 flood: download, channel geometry,
BARC-SS, a rating curve before and after, HAND-FIM, and the flood map. Run it from an
empty folder, since FIMserv writes `code/`, `data/`, and `output/` folders where it runs.

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
- The adjustment is applied at every stage, including stage 0. In BARC (no shift),
  each reach therefore has two stage-0 rows: the zero-flow anchor and the adjusted
  stage-0 row.
- As in the OWP bathymetry adjustment, the subdivision columns (`subdiv_applied`,
  `overbank_n`, `channel_n`, `subdiv_discharge_cms`) are removed from the hydroTable,
  and the original discharge is kept as `discharge_cms_org`.

## Tests

```bash
pip install -e ".[test]"
python -m pytest tests
```

The tests check Equations 1 to 6, that BARC and BARC-SS differ only by the stage
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

## Related repositories

| Repository | What it does |
|------------|--------------|
| [HydroGeom](https://github.com/Reizrb/HydroGeom) | Channel geometry for every NHDPlusV2.1 reach: web page, API, and Python package |
| [Bankfull-and-mean-flow-channel-geometry-for-CONUS](https://github.com/Reizrb/Bankfull-and-mean-flow-channel-geometry-for-CONUS) | The ML models and dataset behind the channel geometry |
| [FIMserv](https://github.com/sdmlua/FIMserv) | Runs OWP HAND-FIM for any HUC8 |

## License

[MIT](LICENSE)

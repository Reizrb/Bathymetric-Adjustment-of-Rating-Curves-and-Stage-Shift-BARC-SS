# BARC-SS: Bathymetric Adjustment of Rating Curves and Stage Shift

DEMs capture the water surface, not the river bed, so the synthetic rating curves of
[OWP HAND-FIM](https://github.com/NOAA-OWP/inundation-mapping) miss the channel below
the water. This code adds that channel to the HAND-FIM hydrotables before the flood map
is made:

- **BARC** adds the bankfull channel area and wetted perimeter to the hydrotable and
  recomputes discharge with Manning's equation.
- **BARC-SS** does the same, then shifts each rating curve's stage up by the reach's
  mean-flow depth and anchors it at zero flow.

Channel width and depth come from [HydroGeom](https://github.com/Reizrb/HydroGeom),
for a channel shape you choose.

<p align="center"><img src="assets/fig4a_step1_initial.png" alt="Initial HAND-FIM rating curve" width="75%"><br><sub>Initial HAND-FIM: the DEM captures only the channel above the water surface.</sub></p>
<p align="center"><img src="assets/fig4b_step2-3_barc.png" alt="BARC" width="75%"><br><sub>BARC: the missing channel geometry is added and discharge recomputed.</sub></p>
<p align="center"><img src="assets/fig4c_step4-5_barc-ss.png" alt="BARC-SS" width="75%"><br><sub>BARC-SS: stage is also shifted by the mean-flow depth and anchored at zero flow.</sub></p>

## Run it with FIMserv_BARC-SS

[FIMserv_BARC-SS](https://github.com/Reizrb/FIMserv_BARC-SS) is a fork of
[FIMserv](https://github.com/sdmlua/FIMserv), which downloads and runs OWP HAND-FIM for
any HUC8, with BARC and BARC-SS added as one extra step. Everything else works the same
as FIMserv. (A pull request to add this step to FIMserv is open.)

```bash
pip install "git+https://github.com/Reizrb/FIMserv_BARC-SS" "git+https://github.com/Reizrb/HydroGeom"
```

```python
import fimserve as fm

huc = "11110203"
fm.DownloadHUC8(huc, version="4.5")
fm.BathymetricAdjustment(huc, method="BARC-SS", shape="a", value=0.4)   # or method="BARC"
fm.getNWMretrospectivedata(huc, "2015-12-26", "2016-01-05", value_times=["2016-01-03 00:00:00"])
fm.runOWPHANDFIM(huc)
```

**Example notebook:** [`examples/BARC-SS_HUC8_11110203_a0.4.ipynb`](examples/BARC-SS_HUC8_11110203_a0.4.ipynb)
runs this for HUC8 11110203 (Arkansas River), a trapezoidal channel with a = 0.4, and the
flood of 3 January 2016, then plots a rating curve before and after and maps the flooded
area. Run it from an empty folder.

| Setting | Options |
|---------|---------|
| `method` | `"BARC"` or `"BARC-SS"` |
| `shape="a"` | trapezoid, `value` = bottom width / top width (0 triangle, 1 rectangle) |
| `shape="r"` | Dingman's r, `value` = r (1 triangle, 2 parabola, larger toward rectangle) |
| `shape="z"` | side slope, `value` = z (z:1) |

The original hydrotables are saved in `output/flood_<huc>/original_hydrotables/` the first
time, and every run starts from them. `fm.RestoreHydrotable(huc)` puts them back, and
`fm.plotBathymetricSRC(huc)` compares a reach's rating curve before and after.

## Use the code on its own

[`barc_ss.py`](barc_ss.py) is the whole method in one file. It works on any HUC8 folder
downloaded with FIMserv:

```python
import barc_ss

flood_dir = "output/flood_11110203"
ids = barc_ss.feature_ids(flood_dir, "11110203")
geometry = barc_ss.channel_geometry(ids, shape="a", value=0.4)
barc_ss.apply_barc(flood_dir, "11110203", geometry, method="BARC-SS")
```

It needs `pandas`, `numpy`, `matplotlib`, and HydroGeomKit
(`pip install "git+https://github.com/Reizrb/HydroGeom"`).

## What it changes

For every stage of each reach's hydrotable, with bankfull area A<sub>bnk</sub> and
wetted perimeter P<sub>bnk</sub> for the chosen shape and reach length L (km):

| | |
|---|---|
| Volume += A<sub>bnk</sub> × L × 1000 | (1) |
| Bed area += P<sub>bnk</sub> × L × 1000 | (2) |
| Wetted perimeter P += P<sub>bnk</sub> | (3) |
| Wet area A += A<sub>bnk</sub> | (4) |
| R = A / P | (5) |
| Q = (1/n) × A × R<sup>2/3</sup> × S<sup>1/2</sup>, with n and S unchanged | (6) |
| BARC-SS only: stage += mean-flow depth | (7) |

Both methods then add a point at stage = 0, discharge = 0 at the start of each rating
curve. The HAND raster is not changed, reaches without a geometry prediction are left
as they are, and the original discharge is kept in `discharge_cms_org`.

## Citation

- Zarrabi, R., Cohen, S., Pruitt, C., Baruah, A., McDermott, R., & Chen, Y. (2026).
  Sensitivity of terrain-based flood inundation model (OWP HAND-FIM) predictions to
  channel geometry: Insights from bathymetric adjustments of rating curves and stage
  shift. *Journal of Hydrology* (under review).
- Zarrabi, R., McDermott, R., Erfani, S. M. H., & Cohen, S. (2025). Bankfull and
  mean-flow channel geometry estimation through machine learning algorithms across
  the CONtiguous United States (CONUS). *Water Resources Research*, 61(2).
  https://doi.org/10.1029/2024WR037997

## License

[MIT](LICENSE)

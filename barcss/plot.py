"""Plots: rating curves before and after adjustment, and the flood map."""
from __future__ import annotations

import glob
import os

import numpy as np
import pandas as pd

from .core import _huc_dir, original_hydrotable


def _applied_method(output_dir, huc8):
    p = os.path.join(output_dir, f"flood_{huc8}", "barcss_applied.txt")
    return open(p).read().strip() if os.path.exists(p) else "Adjusted"


def rating_curve(output_dir, huc8, hydroid=None, ax=None):
    """Plot the original and adjusted rating curve of one reach (HydroID).

    With hydroid=None, picks the reach with the largest mean-flow depth.
    """
    import matplotlib.pyplot as plt
    huc8 = str(huc8)
    ht = pd.read_csv(os.path.join(_huc_dir(output_dir, huc8), "hydrotable.csv"))
    if hydroid is None:
        hydroid = ht.sort_values("shift_datum_m", na_position="first").iloc[-1]["HydroID"]
    new = ht[ht["HydroID"] == hydroid]
    org = original_hydrotable(output_dir, huc8, new["branch_id"].iloc[0])
    org = org[org["HydroID"] == hydroid]

    if ax is None:
        _, ax = plt.subplots(figsize=(6, 4.5))
    ax.plot(org["discharge_cms"], org["stage"], "-", color="#c0504d", lw=2, label="Original HAND-FIM")
    ax.plot(new["discharge_cms"], new["stage"], "-", color="#2f6f9f", lw=2,
            label=_applied_method(output_dir, huc8))
    ax.set_xlabel("Discharge (m³/s)")
    ax.set_ylabel("Stage (m)")
    ax.set_title(f"HydroID {int(hydroid)} (feature_id {int(new['feature_id'].iloc[0])})")
    ax.grid(alpha=.3)
    ax.legend()
    return ax


def inundation_file(output_dir, huc8):
    """Path of the most recent FIMserv inundation raster for the HUC8."""
    files = glob.glob(os.path.join(output_dir, f"flood_{huc8}", f"{huc8}_inundation", "*_inundation.tif"))
    if not files:
        raise FileNotFoundError("No inundation map found. Run fimserve.runOWPHANDFIM(huc) first.")
    return max(files, key=os.path.getmtime)


def flood_map(output_dir, huc8, tif=None, ax=None, color="#1f78b4", title=None):
    """Map only the inundated cells, with the HUC8 boundary for reference."""
    import matplotlib.pyplot as plt
    import rasterio
    from matplotlib.colors import ListedColormap

    huc8 = str(huc8)
    tif = tif or inundation_file(output_dir, huc8)
    with rasterio.open(tif) as r:
        data = r.read(1)
        b, crs, res = r.bounds, r.crs, r.res
    wet = data > 0
    if ax is None:
        _, ax = plt.subplots(figsize=(9, 7))
    ax.imshow(np.ma.masked_where(~wet, np.ones(data.shape)), cmap=ListedColormap([color]),
              extent=(b.left, b.right, b.bottom, b.top), interpolation="nearest", zorder=2)

    boundary = os.path.join(_huc_dir(output_dir, huc8), "wbd8_clp.gpkg")
    if os.path.exists(boundary):
        import geopandas as gpd
        gpd.read_file(boundary).to_crs(crs).boundary.plot(ax=ax, color="0.35", lw=.8, zorder=3)

    area_km2 = wet.sum() * abs(res[0] * res[1]) / 1e6
    ax.set_title(title or f"{_applied_method(output_dir, huc8)} flood map, HUC8 {huc8} "
                          f"({area_km2:,.1f} km² inundated)")
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_aspect("equal")
    return ax

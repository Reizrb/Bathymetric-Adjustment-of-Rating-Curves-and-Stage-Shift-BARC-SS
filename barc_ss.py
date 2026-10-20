"""
BARC and BARC-SS: bathymetric adjustment of OWP HAND-FIM rating curves.

    BARC     adds the bankfull channel area and wetted perimeter missing from the DEM
             to the hydrotable and recomputes discharge with Manning's equation
    BARC-SS  BARC, plus shifts each rating curve's stage up by the mean-flow depth
             and re-anchors it at zero flow

Use on any HUC8 downloaded with FIMserv (https://github.com/sdmlua/FIMserv):

    import barc_ss
    flood_dir = "output/flood_11110203"
    ids = barc_ss.feature_ids(flood_dir, "11110203")
    geometry = barc_ss.channel_geometry(ids, shape="a", value=0.4)
    barc_ss.apply_barc(flood_dir, "11110203", geometry, method="BARC-SS")

Requires pandas, numpy, matplotlib, and HydroGeomKit
(pip install "git+https://github.com/Reizrb/HydroGeom"); the flood map also needs
rasterio, geopandas, and contextily.
"""

import os
import shutil

import numpy as np
import pandas as pd

METHODS = ("BARC", "BARC-SS")
SHAPES = {
    "a": "trapezoid, bottom/top width ratio",
    "r": "Dingman's r",
    "z": "side slope z:1",
}


def channel_geometry(feature_ids, shape="a", value=0.4):
    """Bankfull geometry for each reach from HydroGeom (Zarrabi et al., 2025).

    Gets bankfull width and depth and mean-flow depth for the given NWM feature IDs
    (NHDPlusV2.1 COMIDs), then computes the bankfull cross-sectional area and wetted
    perimeter for the channel shape:
        shape="a": trapezoid, value = bottom width / top width (0 triangle, 1 rectangle)
        shape="r": Dingman's r, value = r (1 triangle, 2 parabola, large = rectangle)
        shape="z": side slope, value = z (z:1)
    """
    try:
        from hydrogeomkit import get_channel_geometry, wetted_perimeter, xsec_area
    except ImportError as e:
        raise ImportError(
            "BARC needs HydroGeomKit for channel geometry. Install it with:\n"
            '    pip install "git+https://github.com/Reizrb/HydroGeom"'
        ) from e
    if shape not in SHAPES:
        raise ValueError(f'shape must be "a", "r", or "z", not {shape!r}')
    geom = get_channel_geometry(comids=list(feature_ids))
    geom = geom.rename(columns={"comid": "feature_id"})
    geom["bnk_xsce_A"] = xsec_area(
        geom["bnk_width"], geom["bnk_depth"], shape=shape, **{shape: value}
    )
    geom["bnk_wet_p"] = wetted_perimeter(
        geom["bnk_width"], geom["bnk_depth"], shape=shape, **{shape: value}
    )
    return geom[
        [
            "feature_id",
            "bnk_width",
            "bnk_depth",
            "mf_width",
            "mf_depth",
            "bnk_xsce_A",
            "bnk_wet_p",
        ]
    ]


def adjust_src(src, geometry):
    """BARC Eqs. 1-6 on one branch's src_full_crosswalked table."""
    src = src.rename(columns={"Discharge (m3s-1)": "Discharge (m3s-1)_Org"})
    src = src.merge(geometry, how="left", on="feature_id")
    src = src.rename(
        columns={
            "bnk_width": "top_width",
            "bnk_xsce_A": "missing_xs_area_m2",
            "bnk_wet_p": "missing_wet_perimeter_m",
            "mf_depth": "shift_datum_m",
        }
    )
    src["Bathymetry_source"] = "Zarrabi_etal"
    src["BARC"] = "TRUE"
    # reaches without a geometry prediction are left unadjusted
    cols = ["missing_xs_area_m2", "missing_wet_perimeter_m", "shift_datum_m"]
    src[cols] = src[cols].fillna(0.0)

    L = src["LENGTHKM"] * 1000
    src["Volume (m3)"] += src["missing_xs_area_m2"] * L  # Eq. 1
    src["BedArea (m2)"] += src["missing_wet_perimeter_m"] * L  # Eq. 2
    src["WettedPerimeter (m)"] += src["missing_wet_perimeter_m"]  # Eq. 3
    src["WetArea (m2)"] += src["missing_xs_area_m2"]  # Eq. 4
    src["HydraulicRadius (m)"] = (
        src["WetArea (m2)"] / src["WettedPerimeter (m)"]
    ).fillna(0)  # Eq. 5
    src["Discharge (m3s-1)_BARC"] = (
        src["WetArea (m2)"]
        * src["HydraulicRadius (m)"] ** (2 / 3)
        * src["SLOPE"] ** 0.5
        / src["ManningN"]
    )  # Eq. 6
    return src


def adjust_hydrotable(ht, src, method="BARC-SS"):
    """Put the BARC discharge in one branch's hydroTable; BARC-SS also shifts stage (Eq. 7)."""
    new = src[
        [
            "HydroID",
            "Stage",
            "Bathymetry_source",
            "Discharge (m3s-1)_BARC",
            "shift_datum_m",
        ]
    ]
    new = new.rename(columns={"Stage": "stage"}).copy()
    new["discharge_cms"] = new["Discharge (m3s-1)_BARC"]

    ht = ht.copy()
    ht["discharge_cms_org"] = ht["discharge_cms"]
    ht = ht.drop(
        columns=[
            "subdiv_applied",
            "discharge_cms",
            "overbank_n",
            "channel_n",
            "subdiv_discharge_cms",
            "Bathymetry_source",
        ],
        errors="ignore",
    )
    ht = ht.merge(new, how="left", on=["HydroID", "stage"])

    if method == "BARC-SS":
        ht["stage"] = ht["stage"] + ht["shift_datum_m"]  # Eq. 7

    # zero-flow anchor: stage = 0, discharge = 0 at the start of every rating curve
    ht["_reach"] = pd.factorize(ht["HydroID"])[0]
    ht["_row"] = np.arange(len(ht))
    anchor = ht.groupby("_reach", sort=False).head(1).copy()
    anchor[["stage", "discharge_cms", "Discharge (m3s-1)_BARC"]] = 0
    anchor["_row"] = -1
    return (
        pd.concat([anchor, ht])
        .sort_values(["_reach", "_row"], kind="stable")
        .drop(columns=["_reach", "_row"])
        .reset_index(drop=True)
    )


def _backup(flood_dir, path):
    """Untouched copy of a HAND-FIM file, made on first use, so reruns start from the original."""
    copy = os.path.join(
        flood_dir, "original_hydrotables", os.path.relpath(path, flood_dir)
    )
    if not os.path.exists(copy):
        os.makedirs(os.path.dirname(copy), exist_ok=True)
        shutil.copy2(path, copy)
    return copy


def feature_ids(flood_dir, huc):
    """NWM feature IDs in the HUC8's rating curves."""
    ids = set()
    for b, src_file, _ in _branches(flood_dir, huc):
        f = os.path.join(
            flood_dir, "original_hydrotables", os.path.relpath(src_file, flood_dir)
        )
        ids.update(
            pd.read_csv(f if os.path.exists(f) else src_file, usecols=["feature_id"])[
                "feature_id"
            ]
            .dropna()
            .astype(int)
        )
    return sorted(ids)


def _branches(flood_dir, huc):
    root = os.path.join(flood_dir, str(huc), "branches")
    if not os.path.isdir(root):
        raise FileNotFoundError(
            f"No HAND-FIM branches at {root}. Download the HUC8 first."
        )
    out = []
    for b in sorted(os.listdir(root)):
        src = os.path.join(root, b, f"src_full_crosswalked_{b}.csv")
        ht = os.path.join(root, b, f"hydroTable_{b}.csv")
        if os.path.exists(src) and os.path.exists(ht):
            out.append((b, src, ht))
    return out


def apply_barc(flood_dir, huc, geometry, method="BARC-SS"):
    """Apply BARC or BARC-SS to every branch of a HUC8 and save the hydrotables in place.

    flood_dir: the HUC8's FIMserv folder, e.g. "output/flood_11110203"
    """
    method = str(method).upper()
    if method not in METHODS:
        raise ValueError(f'method must be "BARC" or "BARC-SS", not {method!r}')
    geometry = geometry.drop_duplicates("feature_id")
    parts = []
    for b, src_file, ht_file in _branches(flood_dir, huc):
        src = adjust_src(pd.read_csv(_backup(flood_dir, src_file)), geometry)
        src.to_csv(src_file, index=False)
        ht = pd.read_csv(
            _backup(flood_dir, ht_file),
            dtype={
                "HUC": str,
                "last_updated": object,
                "submitter": object,
                "obs_source": object,
            },
        )
        ht = adjust_hydrotable(ht, src, method)
        ht.to_csv(ht_file, index=False)
        parts.append(ht.assign(branch_id=b))
    hydrotable = pd.concat(parts, ignore_index=True)
    cols = list(hydrotable.columns)
    cols.insert(1, cols.pop(cols.index("branch_id")))
    top = os.path.join(flood_dir, str(huc), "hydrotable.csv")
    if os.path.exists(top):
        _backup(flood_dir, top)
    hydrotable[cols].to_csv(top, index=False)
    print(f"{method} applied to {len(parts)} branches of HUC8 {huc}.")
    return hydrotable[cols]


def restore_hydrotables(flood_dir):
    """Put back the original HAND-FIM hydrotables."""
    root = os.path.join(flood_dir, "original_hydrotables")
    if not os.path.isdir(root):
        print("Nothing to restore.")
        return
    for d, _, files in os.walk(root):
        for f in files:
            shutil.copy2(
                os.path.join(d, f),
                os.path.join(flood_dir, os.path.relpath(os.path.join(d, f), root)),
            )
    print("Original hydrotables restored.")


def plot_rating_curve(
    flood_dir, huc, hydroid, label="BARC-SS", discharge_csv=None, ax=None
):
    """Original (NO_BARC) and adjusted rating curve of one reach.

    discharge_csv: a FIMserv discharge file (feature_id, discharge); if given, the
    reach's flood discharge is drawn as a dotted vertical line.
    """
    import matplotlib.pyplot as plt

    ht = pd.read_csv(os.path.join(flood_dir, str(huc), "hydrotable.csv"))
    new = ht[ht["HydroID"] == hydroid]
    if new.empty:
        raise ValueError(f"HydroID {hydroid} is not in HUC8 {huc}.")
    b, fid = str(new["branch_id"].iloc[0]), int(new["feature_id"].iloc[0])
    org = pd.read_csv(
        os.path.join(
            flood_dir,
            "original_hydrotables",
            str(huc),
            "branches",
            b,
            f"hydroTable_{b}.csv",
        )
    )
    org = org[org["HydroID"] == hydroid]

    if ax is None:
        _, ax = plt.subplots(figsize=(7, 6))
    ax.set_facecolor("#f7f7f7")
    ax.plot(org["discharge_cms"], org["stage"], color="black", lw=1.8, label="NO_BARC")
    ax.plot(new["discharge_cms"], new["stage"], color="#3a86ff", lw=2.2, label=label)
    if discharge_csv:
        q = pd.read_csv(discharge_csv)
        q = q.loc[q["feature_id"].astype(int) == fid, "discharge"]
        if len(q):
            ax.axvline(
                float(q.iloc[0]), color="black", ls=":", lw=2, label="Flood Discharge"
            )
    ax.set_xlabel("Discharge (cms)", fontsize=14)
    ax.set_ylabel("Stage (m)", fontsize=14)
    ax.set_title(
        f"Synthetic Rating Curve for HydroID={int(hydroid)}\nin HUC8={huc}",
        fontsize=15,
        fontweight="bold",
    )
    ax.ticklabel_format(axis="x", style="sci", scilimits=(0, 0), useMathText=True)
    ax.tick_params(labelsize=12)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.grid(color="0.85", lw=0.8)
    ax.legend(loc="lower right", frameon=False, fontsize=12)
    return ax


def _read_flooded(path, max_px=2500):
    """Flooded cells (value > 0) of a FIM raster, reduced to at most max_px wide.
    A coarse cell counts as flooded if any cell inside it is, so narrow rivers stay visible."""
    import rasterio
    from rasterio.windows import Window

    with rasterio.open(path) as r:
        s = max(1, -(-max(r.width, r.height) // max_px))
        h, w = r.height // s, r.width // s
        wet = np.zeros((h, w), bool)
        step = 256 * s
        for row in range(0, h * s, step):
            n = min(step, h * s - row)
            blk = r.read(1, window=Window(0, row, w * s, n)) > 0
            wet[row // s : (row + n) // s] = blk.reshape(n // s, s, w, s).any(
                axis=(1, 3)
            )
        b = r.bounds
        extent = (b.left, b.left + w * s * r.res[0], b.top - h * s * r.res[1], b.top)
        return wet, extent, r.crs


BASEMAPS = {
    "A": "Esri.WorldImagery",  # satellite
    "B": "Esri.WorldTopoMap",  # topographic
    "C": "CartoDB.Positron",  # light gray
    "D": "CartoDB.Voyager",  # soft street map
    "E": "USGS.USImageryTopo",  # USGS satellite with topo labels
    "F": "Esri.WorldShadedRelief",  # terrain hillshade
}


def plot_flood_map(
    fim_tif,
    boundary_file,
    huc=None,
    basemap="Esri.WorldImagery",
    color="#00a8ff",
    alpha=0.35,
    title=None,
    ax=None,
):
    """Flooded area (semi-transparent) and HUC8 boundary on a web basemap (needs contextily).

    boundary_file: the HUC8 polygon, e.g. output/flood_<huc>/<huc>/wbd.gpkg from FIMserv
    (not wbd8_clp.gpkg, which is the buffered processing extent).
    """
    import contextily as ctx
    import geopandas as gpd
    import matplotlib.pyplot as plt
    import xyzservices.providers as xyz
    from matplotlib.colors import to_rgba
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    wet, extent, crs = _read_flooded(fim_tif)
    rgba = np.zeros(wet.shape + (4,))
    rgba[wet] = to_rgba(color, alpha)
    boundary = gpd.read_file(boundary_file)
    if huc is not None:  # keep only this HUC8 if the file holds several
        for col in ("HUC8", "huc8"):
            if col in boundary.columns:
                boundary = boundary[boundary[col].astype(str) == str(huc)]
    boundary = boundary.dissolve().to_crs(crs)
    line = "white" if "Imagery" in basemap else "black"

    if ax is None:
        _, ax = plt.subplots(figsize=(10, 9))
    ax.imshow(rgba, extent=extent, interpolation="nearest", zorder=3)
    boundary.boundary.plot(ax=ax, color=line, lw=1.6, zorder=4)
    x0, y0, x1, y1 = boundary.total_bounds
    pad = 0.03 * max(x1 - x0, y1 - y0)
    ax.set_xlim(x0 - pad, x1 + pad)
    ax.set_ylim(y0 - pad, y1 + pad)
    ctx.add_basemap(
        ax, crs=crs, source=xyz.query_name(basemap), attribution_size=6, zorder=0
    )
    ax.legend(
        handles=[
            Patch(
                facecolor=to_rgba(color, alpha), edgecolor=color, label="Flooded area"
            ),
            Line2D([], [], color=line, lw=1.6, label="HUC8 boundary"),
        ],
        loc="lower left",
        fontsize=11,
        framealpha=0.85,
    )
    if title:
        ax.set_title(title, fontsize=14, fontweight="bold")
    ax.set_axis_off()
    return ax


class quiet:
    """Hide all printed output, including from tools FIMserv runs in the background
    (for example the AWS download progress). Errors are still raised."""

    def __enter__(self):
        import contextlib
        import logging
        import sys
        import warnings

        self._stack = contextlib.ExitStack()
        devnull = self._stack.enter_context(open(os.devnull, "w"))
        for stream in (sys.stdout, sys.stderr):
            with contextlib.suppress(AttributeError, ValueError, OSError):
                stream.flush()
        self._saved = []
        for fd in (1, 2):
            try:
                self._saved.append((fd, os.dup(fd)))
                os.dup2(devnull.fileno(), fd)
            except OSError:
                pass
        self._stack.enter_context(contextlib.redirect_stdout(devnull))
        self._stack.enter_context(contextlib.redirect_stderr(devnull))
        self._stack.enter_context(warnings.catch_warnings())
        warnings.simplefilter("ignore")
        self._log = logging.root.manager.disable
        logging.disable(logging.CRITICAL)
        return self

    def __exit__(self, *exc):
        import logging

        logging.disable(self._log)
        for fd, saved in self._saved:
            os.dup2(saved, fd)
            os.close(saved)
        self._stack.close()
        return False

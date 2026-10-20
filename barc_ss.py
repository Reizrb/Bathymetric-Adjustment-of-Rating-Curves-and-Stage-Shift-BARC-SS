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
(pip install "git+https://github.com/Reizrb/HydroGeom").
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


def plot_rating_curve(flood_dir, huc, hydroid=None, label="Adjusted"):
    """Original vs adjusted rating curve of one reach (default: deepest mean-flow depth)."""
    import matplotlib.pyplot as plt

    ht = pd.read_csv(os.path.join(flood_dir, str(huc), "hydrotable.csv"))
    if hydroid is None:
        hydroid = ht.loc[ht["shift_datum_m"].idxmax(), "HydroID"]
    new = ht[ht["HydroID"] == hydroid]
    b = str(new["branch_id"].iloc[0])
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
    _, ax = plt.subplots(figsize=(6, 4.5))
    ax.plot(
        org["discharge_cms"],
        org["stage"],
        color="#c0504d",
        lw=2,
        label="Original HAND-FIM",
    )
    ax.plot(new["discharge_cms"], new["stage"], color="#2f6f9f", lw=2, label=label)
    ax.set(
        xlabel="Discharge (m³/s)",
        ylabel="Stage (m)",
        title=f"Rating curve, HydroID {int(hydroid)}",
    )
    ax.grid(alpha=0.3)
    ax.legend()
    return ax

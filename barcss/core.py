"""BARC and BARC-SS: bathymetric adjustment of OWP HAND-FIM rating curves.

Both methods add the channel geometry that the DEM cannot see (bankfull
cross-sectional area and wetted perimeter) to each reach's synthetic rating
curve and recompute discharge with Manning's equation. BARC-SS also shifts the
stage datum up by the mean-flow depth so stage is measured from the channel bed.

    BARC     geometry adjustment only
    BARC-SS  geometry adjustment + stage shift by mean-flow depth
"""
from __future__ import annotations

import os
import shutil
from typing import Union

import numpy as np
import pandas as pd

METHODS = ("BARC", "BARC-SS")
BACKUP_DIR = "barcss_original"   # untouched HAND-FIM files, kept so runs can be repeated


# --------------------------------------------------------------------- paths
def _huc_dir(output_dir: str, huc8: str) -> str:
    return os.path.join(output_dir, f"flood_{huc8}", str(huc8))


def _branches(output_dir: str, huc8: str) -> list:
    path = os.path.join(_huc_dir(output_dir, huc8), "branches")
    if not os.path.isdir(path):
        raise FileNotFoundError(
            f"No HAND-FIM branches found at {path}. Download the HUC8 first, "
            f"for example with fimserve.DownloadHUC8('{huc8}', version='4.5').")
    return sorted(f for f in os.listdir(path) if os.path.isdir(os.path.join(path, f)))


def _branch_files(output_dir, huc8, branch):
    d = os.path.join(_huc_dir(output_dir, huc8), "branches", branch)
    return (os.path.join(d, f"src_full_crosswalked_{branch}.csv"),
            os.path.join(d, f"hydroTable_{branch}.csv"))


def _backup_path(output_dir, huc8, path):
    rel = os.path.relpath(path, _huc_dir(output_dir, huc8))
    return os.path.join(output_dir, f"flood_{huc8}", BACKUP_DIR, rel)


def _original(output_dir, huc8, path):
    """Path of the untouched copy of `path`, creating the copy on first use."""
    backup = _backup_path(output_dir, huc8, path)
    if not os.path.exists(backup):
        os.makedirs(os.path.dirname(backup), exist_ok=True)
        shutil.copy2(path, backup)
    return backup


# ------------------------------------------------------------------ geometry
def feature_ids(output_dir: str, huc8: Union[str, int]) -> list:
    """NWM feature IDs (NHDPlusV2.1 COMIDs) used by the HUC8's rating curves."""
    huc8 = str(huc8)
    ids = set()
    for b in _branches(output_dir, huc8):
        src, _ = _branch_files(output_dir, huc8, b)
        src = _backup_path(output_dir, huc8, src) if os.path.exists(_backup_path(output_dir, huc8, src)) else src
        if os.path.exists(src):
            ids.update(pd.read_csv(src, usecols=["feature_id"])["feature_id"].dropna().astype(int))
    return sorted(ids)


def bankfull_geometry(comids, shape: str = "a", **shape_param) -> pd.DataFrame:
    """Bankfull and mean-flow geometry for each reach, from HydroGeom.

    Gets ML-predicted bankfull and mean-flow width and depth (Zarrabi et al., 2025)
    with HydroGeomKit, then computes the bankfull cross-sectional area and wetted
    perimeter for the chosen channel shape.

    shape: "a" (trapezoid, bottom/top width ratio a), "r" (Dingman's r), or "z" (side slope z:1)
    shape_param: the value for that shape, for example a=0.4

    Returns columns: feature_id, bnk_width, bnk_depth, mf_width, mf_depth,
    bnk_xsce_A (m2), bnk_wet_p (m).
    """
    from hydrogeomkit import get_channel_geometry, xsec_area, wetted_perimeter

    df = get_channel_geometry(comids=list(comids))
    df = df.rename(columns={"comid": "feature_id"})
    df["bnk_xsce_A"] = xsec_area(df["bnk_width"], df["bnk_depth"], shape=shape, **shape_param)
    df["bnk_wet_p"] = wetted_perimeter(df["bnk_width"], df["bnk_depth"], shape=shape, **shape_param)
    cols = ["feature_id", "bnk_width", "bnk_depth", "mf_width", "mf_depth", "bnk_xsce_A", "bnk_wet_p"]
    return df[cols]


# ----------------------------------------------------------------- algorithm
def adjust_src(src_df: pd.DataFrame, geometry: pd.DataFrame) -> pd.DataFrame:
    """Add the missing channel geometry to one branch's SRC table (Eqs. 1-6)."""
    src = src_df.rename(columns={"Discharge (m3s-1)": "Discharge (m3s-1)_Org"})
    src = src.merge(geometry, how="left", on="feature_id")
    src = src.rename(columns={
        "bnk_width": "top_width",
        "bnk_xsce_A": "missing_xs_area_m2",
        "bnk_wet_p": "missing_wet_perimeter_m",
        "mf_depth": "shift_datum_m",
    })
    src["Bathymetry_source"] = "Zarrabi_etal"
    src["BARC"] = "TRUE"
    miss = ["missing_xs_area_m2", "missing_wet_perimeter_m"]
    src[miss] = src[miss].fillna(0.0)
    # Reaches without a geometry prediction are left unadjusted: no added area, no shift.
    src["shift_datum_m"] = src["shift_datum_m"].fillna(0.0)

    length_m = src["LENGTHKM"] * 1000
    src["Volume (m3)"] = src["Volume (m3)"] + src["missing_xs_area_m2"] * length_m            # Eq. 1
    src["BedArea (m2)"] = src["BedArea (m2)"] + src["missing_wet_perimeter_m"] * length_m      # Eq. 2
    src["WettedPerimeter (m)"] = src["WettedPerimeter (m)"] + src["missing_wet_perimeter_m"]   # Eq. 3
    src["WetArea (m2)"] = src["WetArea (m2)"] + src["missing_xs_area_m2"]                      # Eq. 4

    src["HydraulicRadius (m)"] = (src["WetArea (m2)"] / src["WettedPerimeter (m)"]).fillna(0)  # Eq. 5
    src["Discharge (m3s-1)_BARC"] = (                                                          # Eq. 6
        src["WetArea (m2)"] * src["HydraulicRadius (m)"] ** (2.0 / 3)
        * src["SLOPE"] ** 0.5 / src["ManningN"])
    return src


def adjust_hydrotable(htable: pd.DataFrame, src: pd.DataFrame, method: str = "BARC-SS") -> pd.DataFrame:
    """Write the adjusted discharge into one branch's hydroTable.

    BARC-SS also shifts stage up by the mean-flow depth. Both methods then add a
    zero-flow anchor (stage = 0, discharge = 0) at the start of every reach's curve.
    """
    method = _check_method(method)
    trim = src[["HydroID", "Stage", "Bathymetry_source", "Discharge (m3s-1)_BARC", "shift_datum_m"]].copy()
    trim = trim.rename(columns={"Stage": "stage"})
    trim["discharge_cms"] = trim["Discharge (m3s-1)_BARC"]

    ht = htable.copy()
    ht["discharge_cms_org"] = ht["discharge_cms"]
    ht = ht.drop(columns=["subdiv_applied", "discharge_cms", "overbank_n", "channel_n",
                          "subdiv_discharge_cms", "Bathymetry_source"], errors="ignore")
    ht = ht.merge(trim, how="left", on=["HydroID", "stage"])

    if method == "BARC-SS":
        ht["stage"] = ht["stage"] + ht["shift_datum_m"]          # the only difference from BARC

    # Zero-flow anchor: a copy of each reach's first row with stage and discharge set to 0.
    ht["_reach"] = pd.factorize(ht["HydroID"])[0]
    ht["_row"] = np.arange(len(ht))
    anchor = ht.groupby("_reach", sort=False).head(1).copy()
    anchor[["stage", "discharge_cms", "Discharge (m3s-1)_BARC"]] = 0
    anchor["_row"] = -1
    out = (pd.concat([anchor, ht])
             .sort_values(["_reach", "_row"], kind="stable")
             .drop(columns=["_reach", "_row"])
             .reset_index(drop=True))
    return out


def _check_method(method: str) -> str:
    m = str(method).upper().replace("_", "-")
    if m not in METHODS:
        raise ValueError(f'method must be "BARC" or "BARC-SS", not {method!r}')
    return m


# -------------------------------------------------------------- run on a HUC
def apply(output_dir: str, huc8: Union[str, int], geometry: pd.DataFrame,
          method: str = "BARC-SS", verbose: bool = True) -> pd.DataFrame:
    """Apply BARC or BARC-SS to every branch of a downloaded HAND-FIM HUC8.

    output_dir: the FIMserv "output" folder (holds flood_<huc8>/)
    geometry:   table from bankfull_geometry() (feature_id, bnk_xsce_A, bnk_wet_p, mf_depth)
    method:     "BARC" or "BARC-SS"

    The first run saves untouched copies of the HAND-FIM files in
    flood_<huc8>/barcss_original/. Every run starts from those copies, so you can
    run again with another method or channel shape. Use restore() to go back.
    Returns the combined adjusted hydrotable.
    """
    huc8 = str(huc8)
    method = _check_method(method)
    need = {"feature_id", "bnk_xsce_A", "bnk_wet_p", "mf_depth"}
    if not need.issubset(geometry.columns):
        raise ValueError(f"geometry needs columns {sorted(need)}")
    geometry = geometry.drop_duplicates("feature_id")

    parts = []
    for branch in _branches(output_dir, huc8):
        src_file, ht_file = _branch_files(output_dir, huc8, branch)
        if not (os.path.exists(src_file) and os.path.exists(ht_file)):
            continue
        src = adjust_src(pd.read_csv(_original(output_dir, huc8, src_file)), geometry)
        src.to_csv(src_file, index=False)

        ht = pd.read_csv(_original(output_dir, huc8, ht_file), dtype={
            "HUC": str, "last_updated": object, "submitter": object, "obs_source": object})
        ht = adjust_hydrotable(ht, src, method)
        ht.to_csv(ht_file, index=False)

        ht = ht.copy()
        ht["branch_id"] = branch
        parts.append(ht)

    if not parts:
        raise FileNotFoundError(f"No branch rating-curve files found for HUC8 {huc8}.")
    hydrotable = pd.concat(parts, ignore_index=True)
    cols = list(hydrotable.columns)
    cols.insert(1, cols.pop(cols.index("branch_id")))
    hydrotable = hydrotable[cols]

    top = os.path.join(_huc_dir(output_dir, huc8), "hydrotable.csv")
    if os.path.exists(top):
        _original(output_dir, huc8, top)
    hydrotable.to_csv(top, index=False)
    with open(os.path.join(output_dir, f"flood_{huc8}", "barcss_applied.txt"), "w") as f:
        f.write(f"{method}\n")
    if verbose:
        print(f"{method} applied to {len(parts)} branches of HUC8 {huc8}.")
    return hydrotable


def restore(output_dir: str, huc8: Union[str, int]) -> None:
    """Put back the original HAND-FIM files (undo apply)."""
    huc8 = str(huc8)
    root = os.path.join(output_dir, f"flood_{huc8}", BACKUP_DIR)
    if not os.path.isdir(root):
        print("Nothing to restore: BARC/BARC-SS has not been applied here.")
        return
    for d, _, files in os.walk(root):
        for f in files:
            backup = os.path.join(d, f)
            target = os.path.join(_huc_dir(output_dir, huc8), os.path.relpath(backup, root))
            shutil.copy2(backup, target)
    marker = os.path.join(output_dir, f"flood_{huc8}", "barcss_applied.txt")
    if os.path.exists(marker):
        os.remove(marker)
    print(f"Original HAND-FIM files restored for HUC8 {huc8}.")


def original_hydrotable(output_dir: str, huc8: Union[str, int], branch: str) -> pd.DataFrame:
    """The unadjusted hydroTable of one branch (for comparison plots)."""
    _, ht_file = _branch_files(output_dir, str(huc8), str(branch))
    backup = _backup_path(output_dir, str(huc8), ht_file)
    return pd.read_csv(backup if os.path.exists(backup) else ht_file)


def keep_flood_map(output_dir: str, huc8: Union[str, int], label: str) -> str:
    """Rename the latest FIMserv flood map so it is kept, e.g. label="BARC-SS_a0.4".

    FIMserv skips a run when a flood map for the same discharge file already
    exists. Renaming the map lets you run FIM again with another method or shape.
    """
    import glob
    huc8 = str(huc8)
    d = os.path.join(output_dir, f"flood_{huc8}", f"{huc8}_inundation")
    files = [f for f in glob.glob(os.path.join(d, "*_inundation.tif"))]
    if not files:
        raise FileNotFoundError("No flood map to keep. Run fimserve.runOWPHANDFIM(huc) first.")
    latest = max(files, key=os.path.getmtime)
    new = latest.replace("_inundation.tif", f"_{label}_inundation.tif")
    os.replace(latest, new)
    print(f"Saved as {os.path.basename(new)}")
    return new

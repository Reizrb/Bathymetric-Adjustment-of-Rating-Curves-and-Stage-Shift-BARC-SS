import os
import shutil
import sys
import warnings

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(__file__))
import barcss
from make_fake_huc import HUC, make

BRANCHES = ["0", "2284000046"]


def _ht(out, b):
    return pd.read_csv(os.path.join(out, f"flood_{HUC}", HUC, "branches", b, f"hydroTable_{b}.csv"))


def _src(out, b):
    return pd.read_csv(os.path.join(out, f"flood_{HUC}", HUC, "branches", b, f"src_full_crosswalked_{b}.csv"))


@pytest.fixture
def huc(tmp_path):
    out = str(tmp_path / "output")
    geom = make(out)
    return out, geom


def test_matches_original_code(tmp_path):
    """BARC-SS gives the same rating curves as the code used in the paper."""
    from reference.original_utils import modify_hydrotable_bathymetric_adjustments as original
    a, b = str(tmp_path / "a"), str(tmp_path / "b")
    geom = make(a, with_missing_reach=False); make(b, with_missing_reach=False)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        original(a, HUC, geom.copy())
    barcss.apply(b, HUC, geom, method="BARC-SS", verbose=False)
    for br in BRANCHES:
        old = _ht(a, br).rename(columns={"shift_datume_m": "shift_datum_m"})
        new = _ht(b, br)
        assert list(old.columns) == list(new.columns)
        anchor = new["stage"].eq(0) & new["discharge_cms"].eq(0)
        cols = [c for c in new.columns if c != "Discharge (m3s-1)_BARC"]
        pd.testing.assert_frame_equal(old[cols], new[cols])
        # only cosmetic change: the anchor row's BARC discharge is now 0, as its discharge_cms
        assert (new.loc[anchor, "Discharge (m3s-1)_BARC"] == 0).all()
        pd.testing.assert_series_equal(old.loc[~anchor, "Discharge (m3s-1)_BARC"],
                                       new.loc[~anchor, "Discharge (m3s-1)_BARC"])
        so = _src(a, br).rename(columns={"shift_datume_m": "shift_datum_m"})
        pd.testing.assert_frame_equal(so, _src(b, br))


def test_barc_differs_only_by_stage_shift(tmp_path):
    a, b = str(tmp_path / "a"), str(tmp_path / "b")
    geom = make(a); make(b)
    barcss.apply(a, HUC, geom, method="BARC", verbose=False)
    barcss.apply(b, HUC, geom, method="BARC-SS", verbose=False)
    for br in BRANCHES:
        x, y = _ht(a, br), _ht(b, br)
        pd.testing.assert_frame_equal(x.drop(columns="stage"), y.drop(columns="stage"))
        real = ~(y["stage"].eq(0) & y["discharge_cms"].eq(0) & y.index.isin(y.groupby("HydroID").head(1).index))
        np.testing.assert_allclose(y.loc[real, "stage"], x.loc[real, "stage"] + x.loc[real, "shift_datum_m"])


def test_equations(huc):
    out, geom = huc
    org = _src(out, "0")
    barcss.apply(out, HUC, geom, verbose=False)
    new = _src(out, "0")
    g = geom.set_index("feature_id").loc[new["feature_id"]]
    L = org["LENGTHKM"] * 1000
    np.testing.assert_allclose(new["WetArea (m2)"], org["WetArea (m2)"] + g["bnk_xsce_A"].values)
    np.testing.assert_allclose(new["WettedPerimeter (m)"], org["WettedPerimeter (m)"] + g["bnk_wet_p"].values)
    np.testing.assert_allclose(new["Volume (m3)"], org["Volume (m3)"] + g["bnk_xsce_A"].values * L)
    np.testing.assert_allclose(new["BedArea (m2)"], org["BedArea (m2)"] + g["bnk_wet_p"].values * L)
    R = new["WetArea (m2)"] / new["WettedPerimeter (m)"]
    Q = new["WetArea (m2)"] * R ** (2 / 3) * new["SLOPE"] ** .5 / new["ManningN"]
    np.testing.assert_allclose(new["Discharge (m3s-1)_BARC"], Q)


def test_rerun_does_not_double_count(huc):
    out, geom = huc
    barcss.apply(out, HUC, geom, verbose=False)
    first = _ht(out, "0")
    barcss.apply(out, HUC, geom, verbose=False)
    pd.testing.assert_frame_equal(first, _ht(out, "0"))


def test_switch_method_and_restore(huc):
    out, geom = huc
    original = _ht(out, "0")
    barcss.apply(out, HUC, geom, method="BARC-SS", verbose=False)
    barcss.apply(out, HUC, geom, method="BARC", verbose=False)
    assert "BARC" == open(os.path.join(out, f"flood_{HUC}", "barcss_applied.txt")).read().strip()
    barcss.restore(out, HUC)
    pd.testing.assert_frame_equal(original, _ht(out, "0"))


def test_reach_without_geometry_left_unchanged(huc):
    out, geom = huc
    org = _ht(out, "2284000046")
    barcss.apply(out, HUC, geom, verbose=False)
    new = _ht(out, "2284000046")
    missing = set(org["feature_id"]) - set(geom["feature_id"])
    assert missing
    rows = new[new["feature_id"].isin(missing)].iloc[1:]      # skip the zero anchor
    o = org[org["feature_id"].isin(missing)]
    np.testing.assert_allclose(rows["stage"], o["stage"])      # no shift, no NaN stages
    assert not new["stage"].isna().any()


def test_feature_ids_and_bad_method(huc):
    out, geom = huc
    assert len(barcss.feature_ids(out, HUC)) == 6
    with pytest.raises(ValueError):
        barcss.apply(out, HUC, geom, method="SS")
    with pytest.raises(FileNotFoundError):
        barcss.apply(out, "99999999", geom)


def test_plots(huc, tmp_path):
    import matplotlib
    matplotlib.use("Agg")
    import rasterio
    from rasterio.transform import from_origin
    out, geom = huc
    barcss.apply(out, HUC, geom, verbose=False)
    barcss.plot.rating_curve(out, HUC)
    d = os.path.join(out, f"flood_{HUC}", f"{HUC}_inundation"); os.makedirs(d)
    arr = np.zeros((50, 50), "uint8"); arr[20:30, :] = 1
    with rasterio.open(os.path.join(d, "NWM_x_inundation.tif"), "w", driver="GTiff", height=50, width=50,
                       count=1, dtype="uint8", crs="EPSG:5070", transform=from_origin(0, 500, 10, 10)) as r:
        r.write(arr, 1)
    ax = barcss.plot.flood_map(out, HUC)
    assert "0.0" in ax.get_title() or "km²" in ax.get_title()


def test_keep_flood_map(huc):
    out, geom = huc
    d = os.path.join(out, f"flood_{HUC}", f"{HUC}_inundation"); os.makedirs(d)
    open(os.path.join(d, f"NWM_20160103000000_{HUC}_inundation.tif"), "w").close()
    new = barcss.keep_flood_map(out, HUC, "BARC-SS_a0.4")
    assert new.endswith(f"NWM_20160103000000_{HUC}_BARC-SS_a0.4_inundation.tif")
    assert not os.path.exists(os.path.join(d, f"NWM_20160103000000_{HUC}_inundation.tif"))

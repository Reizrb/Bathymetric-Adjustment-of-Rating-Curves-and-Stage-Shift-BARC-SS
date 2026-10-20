"""Build a small fake HAND-FIM HUC8 folder with the same file layout and columns
as FIMserv downloads, for testing without the real data."""
import os
import numpy as np
import pandas as pd

HUC = "11110203"


def make(output_dir, with_missing_reach=True):
    rng = np.random.default_rng(0)
    geom_rows = []
    hid = 23650001
    for branch in ["0", "2284000046"]:
        bdir = os.path.join(output_dir, f"flood_{HUC}", HUC, "branches", branch)
        os.makedirs(bdir, exist_ok=True)
        src_rows, ht_rows = [], []
        for k in range(3):
            fid = 7787400 + hid % 1000
            length, slope, n = rng.uniform(.3, 2), rng.uniform(1e-4, 3e-3), .06
            for stage in np.round(np.arange(0, 3.05, .3048), 4):
                cells = 0 if stage == 0 else int(200 + 900 * stage)
                top_w = cells * 10 / (length * 1000) if cells else 0.0
                area = top_w * stage * .6
                perim = top_w + 2 * stage
                q = 0.0 if area == 0 else area * (area / perim) ** (2 / 3) * slope ** .5 / n
                src_rows.append({"HydroID": hid, "branch_id": branch, "feature_id": fid, "NextDownID": -1,
                                 "order_": 4, "Number of Cells": cells, "SurfaceArea (m2)": cells * 100.0,
                                 "BedArea (m2)": perim * length * 1000, "TopWidth (m)": top_w,
                                 "LENGTHKM": length, "AREASQKM": 50.0, "WettedPerimeter (m)": perim,
                                 "HydraulicRadius (m)": area / perim if perim else 0.0, "WetArea (m2)": area,
                                 "Volume (m3)": area * length * 1000, "SLOPE": slope, "ManningN": n,
                                 "Stage": stage, "Discharge (m3s-1)": q})
                ht_rows.append({"HydroID": hid, "feature_id": fid, "HUC": HUC, "LakeID": -999,
                                "last_updated": None, "submitter": None, "obs_source": None,
                                "precalb_discharge_cms": q, "calb_coef_final": 1.0, "stage": stage,
                                "discharge_cms": q, "subdiv_applied": False, "overbank_n": .12,
                                "channel_n": .06, "subdiv_discharge_cms": q})
            if not (with_missing_reach and branch == "2284000046" and k == 2):
                w = rng.uniform(30, 90); d = rng.uniform(1.5, 3); mfd = d * rng.uniform(.4, .7)
                geom_rows.append({"feature_id": fid, "bnk_width": w, "bnk_depth": d,
                                  "mf_width": w * .9, "mf_depth": mfd})
            hid += 1
        pd.DataFrame(src_rows).to_csv(os.path.join(bdir, f"src_full_crosswalked_{branch}.csv"), index=False)
        pd.DataFrame(ht_rows).to_csv(os.path.join(bdir, f"hydroTable_{branch}.csv"), index=False)
    pd.DataFrame({"x": [1]}).to_csv(os.path.join(output_dir, f"flood_{HUC}", HUC, "hydrotable.csv"), index=False)
    geom = pd.DataFrame(geom_rows)
    from hydrogeomkit import xsec_area, wetted_perimeter
    geom["bnk_xsce_A"] = xsec_area(geom["bnk_width"], geom["bnk_depth"], shape="a", a=0.4)
    geom["bnk_wet_p"] = wetted_perimeter(geom["bnk_width"], geom["bnk_depth"], shape="a", a=0.4)
    return geom

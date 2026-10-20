"""BARC and BARC-SS: bathymetric adjustment of OWP HAND-FIM rating curves.

    >>> import barcss
    >>> ids = barcss.feature_ids("output", "11110203")
    >>> geom = barcss.bankfull_geometry(ids, shape="a", a=0.4)
    >>> barcss.apply("output", "11110203", geom, method="BARC-SS")
"""
from .core import (METHODS, adjust_hydrotable, adjust_src, apply, bankfull_geometry,
                   feature_ids, keep_flood_map, original_hydrotable, restore)
from . import plot

__version__ = "1.0.0"
__all__ = ["METHODS", "adjust_hydrotable", "adjust_src", "apply", "bankfull_geometry",
           "feature_ids", "keep_flood_map", "original_hydrotable", "restore", "plot"]

#!/usr/bin/env python

# SPDX-FileCopyrightText: 2020-2021 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later

# coding: utf-8

import h5py
import numpy as np

# import tables # Need this, otherwise errors from h5py show up (no idea why; SO be praised)
from pathlib import Path
import pandas as pd
import xarray as xr

from _helpers import configure_logging
import logging

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    configure_logging(snakemake)

    region = snakemake.config["regions"][snakemake.wildcards["region"]]
    year = int(snakemake.wildcards["year"])

    # Different classes (site qualities) of cells for renewables
    r_classes = {
        res: snakemake.config["GlobalEnergyGIS"][res]["classes"]["min"]
        for res in {"wind", "solar"}
    }

    # Dimension names associated with the techs and variables
    dimensions = {
        # wind
        "CFtime_windonshoreA": ["class", "region", "time"],
        "CFtime_windonshoreB": ["class", "region", "time"],
        "CFtime_windoffshore": ["class", "region", "time"],
        "capacity_offshore": ["class", "region"],
        "capacity_onshoreA": ["class", "region"],
        "capacity_onshoreB": ["class", "region"],
        # solar
        "CFtime_cspplantA": ["class", "region", "time"],
        "CFtime_cspplantB": ["class", "region", "time"],
        "CFtime_pvplantA": ["class", "region", "time"],
        "CFtime_pvplantB": ["class", "region", "time"],
        "CFtime_pvrooftop": ["class", "region", "time"],
        "capacity_cspplantA": ["class", "region"],
        "capacity_cspplantB": ["class", "region"],
        "capacity_pvplantA": ["class", "region"],
        "capacity_pvplantB": ["class", "region"],
        "capacity_pvrooftop": ["class", "region"],
        "csp_density": ["dim_0"],
        "pv_density": ["dim_0"],
        "solar_overlap_areaA": ["class_1", "class_2", "region"],
        "solar_overlap_areaB": ["class_1", "class_2", "region"],
    }

    # List as cotainer for single DataArrays
    das = []

    # Open single files consecutively
    for tech in ["wind", "solar"]:
        fn = snakemake.input[tech]
        logger.info(f"Processing: {fn}")

        # Coordinates associated with the different dimensions
        # Classes are different for each tech
        coordinates = {
            "region": sorted(
                list(snakemake.config["regions"][snakemake.wildcards["region"]])
            ),
            "time": pd.date_range(
                start=str(year), end=str(year + 1), freq="H", closed="left"
            ),
            "dim_0": [0],
            "class": r_classes[tech],
            "class_1": r_classes[tech],
            "class_2": r_classes[tech],
        }

        with h5py.File(fn, "r") as f:

            # Read each variable seperatley, convert to DataArray

            for k, v in f.items():
                das.append(
                    xr.DataArray(
                        np.array(v),
                        name=k,
                        dims=dimensions[k],
                        coords={d: coordinates[d] for d in dimensions[k]},
                    )
                )
                logger.info(f"Variable: {k}.")

    # Merge single DataArrays and save Dataset
    logger.info(f"Merging into single dataset and saving.")
    ds = xr.merge(das)
    ds.to_netcdf(snakemake.output[0])

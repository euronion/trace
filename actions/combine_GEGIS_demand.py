#!/usr/bin/env python

# SPDX-FileCopyrightText: 2020-2021 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later

# coding: utf-8

from _helpers import configure_logging

import h5py
import tables  # Required for some mysterious reason for h5py to work properly
import numpy as np
from pathlib import Path
import pandas as pd
import xarray as xr

from _helpers import configure_logging
import logging

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    configure_logging(snakemake)

    region = sorted(
        list(snakemake.config["regions"][snakemake.wildcards["region"]].keys())
    )
    year = int(snakemake.wildcards["era_year"])

    with h5py.File(snakemake.input[0], "r") as f:
        da = xr.DataArray(
            np.array(f["demand"]),
            name="demand",
            dims=["region", "time"],
            coords={
                "region": region,
                "time": pd.date_range(
                    start=str(year), end=str(year + 1), freq="H", closed="left"
                ),
            },
        )

        da.to_netcdf(snakemake.output["netcdf"])
    da = da.sum(dim="time").rename("demand [GWh]") / 1e3
    df = da.to_dataframe()
    df.to_csv(snakemake.output["csv"], sep=";", quotechar='"')

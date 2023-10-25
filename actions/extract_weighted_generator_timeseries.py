#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from pathlib import Path

import hvplot.pandas
import hvplot.xarray
import pypsa
import xarray as xr
from _helpers import configure_logging

if __name__ == "__main__":
    scenario = snakemake.params["scenario"]

    configure_logging(snakemake)

    network = pypsa.Network(snakemake.input.network)

    # Total optimal capacity per carrier
    p_nom_opt = network.generators.groupby("carrier").sum()["p_nom_opt"].to_xarray()
    p_nom_opt = p_nom_opt.rename({"carrier": "technology"})

    total_capacity_per_carrier = network.generators.groupby(
        by="carrier", axis="index"
    ).sum()["p_nom_opt"]

    # Create groupnames without numbers per column to groupby, e.g. "offwind 1" -> "offwind"
    groups = network.generators_t["p_max_pu"].columns.str.replace(
        "\s\d\d?", "", regex=True
    )

    weighted_p_max = (
        (network.generators_t["p_max_pu"] * network.generators["p_nom_opt"])
        .groupby(by=groups, axis="columns")
        .sum()
    )

    weighted_p_max_pu = weighted_p_max / total_capacity_per_carrier

    # Convert structure to nice pandas DataFrame and then to xarray DataArray
    weighted_p_max_pu = (
        weighted_p_max_pu.reset_index()
        .melt("snapshot", var_name="technology", value_name="p_max_pu")
        .set_index(["technology", "snapshot"])
        .to_xarray()
    )

    ds = xr.merge([p_nom_opt, weighted_p_max_pu])

    ds.to_netcdf(snakemake.output[0])


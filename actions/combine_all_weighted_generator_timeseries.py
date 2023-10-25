#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import logging
from pathlib import Path

import xarray as xr
from _helpers import configure_logging

logger = logging.getLogger(__name__)


# In[ ]:


if __name__ == "__main__":
    configure_logging(snakemake)

    ds = []  # Container, used to merge all network-individual results later
    for fp in snakemake.input["networks"]:
        fp = Path(fp)
        logger.info(f"Reading {fp}.")
        ds_network = xr.open_dataset(fp)

        ds_network = ds_network.expand_dims(
            {
                # Add exporter and importer code as new dimension
                "exporter": [
                    "-".join(fp.parts[-2].split("-")[:-1])
                ],  # Allow for "-" in exporter, needs to be re-added
                "importer": [fp.parts[-2].split("-")[-1]],
            }
        )

        # Add individual time-series to list, combined later
        ds.append(ds_network)

    # Combine DataSets
    ds = xr.merge(ds)
    logger.info("Saving combined results...")

    ds.to_netcdf(snakemake.output["combined_timeseries"])


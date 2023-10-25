#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import logging
import pickle
import re

import dateutil
import numpy as np
import pandas as pd
import pypsa
from _helpers import configure_logging

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    scenario = snakemake.params["scenario"]
    configure_logging(snakemake)

    with open(snakemake.input["additional_components"], "rb") as f:
        override_component_attrs = pickle.load(f)

    network = pypsa.Network(override_component_attrs=override_component_attrs)
    network.import_from_netcdf(snakemake.input["network"])

    # Determine all loads related to imports (on the importer side of the ESC)
    import_loads = network.loads[network.loads.bus.str.contains("imp")]

    # Read load profile from file
    df = pd.read_csv(snakemake.input["import_profile"], index_col="snapshot")

    # Conversion to datetime required for pd.update() later
    # Index must match that of network.snapshots
    df.index = pd.to_datetime(df.index)

    # File format consistency checks
    if {"snapshot", "profile [p.u.]"} != set(
        pd.read_csv(snakemake.input["import_profile"]).columns
    ):
        logger.error(
            f"Import profile file '{snakemake.input['import_profile']}' must have "
            f"exactly two columns 'snapshot' and 'profile [p.u.]'"
        )
    if df.index.equals(network.snapshots) is False:
        logger.error(
            f"The 'snapshot' index from '{snakemake.input['import_profile']}' "
            f"does not match the network snapshots: {network.snapshots[:5], ...}. "
            f"One or more wrong entries in the input file?"
        )

    # Create import_profile with each demand attached to an import bus the network
    import_profiles = pd.DataFrame(
        {n: p["p_set"] * df["profile [p.u.]"] for n, p in import_loads.iterrows()}
    )

    # Add import profiles to network
    network.loads_t["p_set"] = network.loads_t["p_set"].join(import_profiles)

    # Remove static demand for these loads
    network.loads.loc[import_profiles.columns, "p_set"] = 0

    # Add optional buffers to the import loads
    if scenario.get("import_buffer", False) is False:
        logger.info("Imports are not buffered.")
    else:
        freqs = {
            "annually": "Y",
            "quaterly": "Q",
            "monthly": "M",
            "weekly": "W",
            "biweekly": "SM",
            "daily": "D",
        }

        if scenario["import_buffer"] not in freqs:
            logger.error(
                f"Unknown 'import_buffer' option '{scenario['import_buffer']}'."
            )

        # Buffer sized to match the maximum demand of any time-group based on the buffer duration,
        # i.e. the buffer size can fully accomodate the demand within the buffer duration
        buffer_e_nom_max = network.loads_t["p_set"][import_loads.index].groupby(
            pd.Grouper(freq=freqs[scenario["import_buffer"]])
        )
        buffer_e_nom_max = np.round(buffer_e_nom_max.sum().max().values)

        logger.info(
            f"Imports are buffered {scenario['import_buffer']}, "
            f"adding {buffer_e_nom_max} MWh buffer to {import_loads.index.values} ."
        )

        # Implement import_buffer by adding a store
        # to the bus of each load on the import side
        # with limited maximum capacity equal to max. time-span ("import_buffer") import demand
        network.madd(
            "Store",
            names="Buffer: " + import_loads.index,
            bus=import_loads.bus.values,
            e_nom_extendable=True,
            e_nom_max=buffer_e_nom_max,
            e_cyclic=True,
        )

    network.export_to_netcdf(snakemake.output["network"])


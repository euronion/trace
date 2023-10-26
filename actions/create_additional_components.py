# SPDX-FileCopyrightText: 2023 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later

import logging
import pickle
from pathlib import Path

import numpy as np
import pypsa
from _helpers import configure_logging

logger = logging.getLogger(__name__)
if __name__ == "__main__":
    configure_logging(snakemake)
    override_component_attrs = pypsa.descriptors.Dict(
        {k: v.copy() for k, v in pypsa.components.component_attrs.items()}
    )

    override_component_attrs["Link"].loc["bus2"] = [
        "string",
        np.nan,
        np.nan,
        "Name of optional 3rd bus to which link is attached.",
        "Input (optional)",
    ]
    override_component_attrs["Link"].loc["efficiency2"] = [
        "static or series",
        "per unit",
        np.nan,
        "Efficiency of power transfer from bus0 to bus2 (static or time-dependent)",
        "Input (optional)",
    ]
    override_component_attrs["Link"].loc["p2"] = [
        "series",
        "MW",
        0.0,
        "3rd bus output",
        "Output",
    ]

    override_component_attrs["Bus"].loc["unit"] = [
        "string",
        np.nan,
        np.nan,
        "Unit of the commoditiy in of the bus. Purely descriptive.",
        "Input (optional)",
    ]

    override_component_attrs["Link"].loc["bus3"] = [
        "string",
        np.nan,
        np.nan,
        "Name of optional 4th bus to which link is attached.",
        "Input (optional)",
    ]
    override_component_attrs["Link"].loc["efficiency3"] = [
        "static or series",
        "per unit",
        np.nan,
        "Efficiency of power transfer from bus0 to bus3 (static or time-dependent)",
        "Input (optional)",
    ]
    override_component_attrs["Link"].loc["p3"] = [
        "series",
        "MW",
        0.0,
        "4th bus output",
        "Output",
    ]

    override_component_attrs["Bus"].loc["unit"] = [
        "string",
        np.nan,
        np.nan,
        "Unit of the commoditiy in of the bus. Purely descriptive.",
        "Input (optional)",
    ]

    output = Path(snakemake.output["additional_components"])

    output.parent.mkdir(parents=True, exist_ok=True)

    with open(snakemake.output["additional_components"], "wb") as f:
        pickle.dump(override_component_attrs, f)

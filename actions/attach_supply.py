# SPDX-FileCopyrightText: 2020-2023 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later

import calendar
import logging
import pickle
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pypsa
import xarray as xr
from _helpers import calculate_annual_investment, configure_logging, read_efficiencies

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    scenario = snakemake.params["scenario"]
    configure_logging(snakemake)

    supply = xr.open_dataset(snakemake.input["supply"])
    logger.info(
        f"Removing the following supply technologies to conform with scenario supply_technologies: {set(supply['technology'].values) - set(scenario['supply_technologies'])}"
    )

    supply = supply.where(lambda x: x.technology.isin(scenario["supply_technologies"]))
    logger.info(f"Remaining supply technologies: {supply['technology'].values}")

    # Load additional components
    with open(snakemake.input["additional_components"], "rb") as f:
        override_component_attrs = pickle.load(f)

    network = pypsa.Network(override_component_attrs=override_component_attrs)
    network.import_from_netcdf(snakemake.input["network"])

    Path(snakemake.output["lcoes"]).parent.mkdir(parents=True, exist_ok=True)

    # Load annual demand (consider overwrite option in config)
    if scenario["synthetic_demand"].lower() in ["custom"]:
        demand = pd.read_csv(
            snakemake.input["demand"], index_col="region", keep_default_na=False
        )

        demand = demand.loc[snakemake.wildcards["from"], "demand [GWh]"]
        demand *= 1e3  # in MWh
        # Apply scenario modifier
        demand *= scenario["modifiers"]["synthetic_demand"]

    else:
        logging.error(
            f"Option {scenario['synthetic_demand']} in config for scenario -> synthethic_demand unknown."
        )

    # Remove capacity classes considered negligible from the supply side.
    # Entries are substituted with "NaN" and then later not considered
    # in the LCoE calculation and subsequently do not enter the model
    # as generators
    _mask = (
        supply["capacities"]
        < snakemake.params["renewable_details"]["ignore_capacities_below"]
    )
    logger.info(
        f"Removing {_mask.sum().item()} renewable classes with individual capacity below "
        f"'ignore_capacity_below'={snakemake.params['renewable_details']['ignore_capacities_below']} MW "
        f"and a combined capacity of {supply['capacities'].where(_mask, 0).sum().item():.2f} MW "
        f"from the available supply."
    )
    supply = supply.where(np.logical_not(_mask))

    capacities = supply["capacities"]
    capacity_factors = supply["profiles"].mean(dim="time")

    # Determine which generation capacities are available for export
    # and which are needed for satisfying domestic electricity demand
    # Simple heuristic:
    # * Consider annual demand and generation per source+source class
    # * Calculate annual costs per technology and LCoE per class/technology
    # * Exclude (fully and partially) sources with lowest LCoE for demand

    wacc = pd.read_csv(
        snakemake.input["wacc"], comment="#", index_col="region", keep_default_na=False
    )
    wacc = wacc.loc[snakemake.wildcards["from"], scenario["wacc"]]
    wacc *= scenario["modifiers"]["wacc"]  # Apply scenario modifier

    def _calculate_investment(name):
        # Auxiliary function required for .map()
        # Calculate annual investment for technology "x"
        # Combine multiple technologies for CSP for estimating LCoE based on exogeneous solar multiple

        if name != "csp-tower":
            invest = calculate_annual_investment(name, wacc, snakemake.input["costs"])
        elif name == "csp-tower":
            # LCoE calculation of CSP needs to account for solar multiple and additional
            # Components for generation of electricity
            solar_multiple = snakemake.params["renewable_details"]["csp_tower"][
                "solar_multiple_estimate"
            ]

            invest = (
                calculate_annual_investment("csp-tower", wacc, snakemake.input["costs"])
                +
                # Share of heat directly converted to electricity
                # Cost in kW_e -> have to consider conversion efficiency for kW_th
                1
                / solar_multiple
                * calculate_annual_investment(
                    "csp-tower power block", wacc, snakemake.input["costs"]
                )
                * read_efficiencies(
                    snakemake.input["efficiencies"], snakemake.wildcards["year"]
                )
                .query("process == 'csp-tower power block'")["efficiency"]
                .item()
                # Store heat not converted directly in TES
                + (solar_multiple - 1)
                / solar_multiple
                * calculate_annual_investment(
                    "csp-tower TES", wacc, snakemake.input["costs"]
                )
            )

        invest *= 1.0e3  # Investment costs in costs data given in kW, here MW is used

        return invest

    ## Annual investment per maximum capacity of each class and technology
    # Investment per MW capacity (use pandas mapping, doesn't seem to work on xarray)
    investment = (
        capacities["technology"].to_pandas().map(_calculate_investment).to_xarray()
    )
    investment = investment * capacities
    investment = investment.rename("annual investment")

    # Calculate annual generation
    generation = (supply["profiles"] * capacities).fillna(0.0).sum(dim="time")
    generation = generation.rename("annual generation")

    # Special treatment if CSP technologies are present
    # Electricity generation from CSP is reduced by efficiency of the CSP plant's power block
    # Reported generation capacity for CSP given in thermal energy, here compare electric energy
    if "csp-tower" in generation.technology:
        generation.loc[dict(technology="csp-tower")] *= (
            read_efficiencies(
                snakemake.input["efficiencies"], snakemake.wildcards["year"]
            )
            .query("process == 'csp-tower power block'")["efficiency"]
            .item()
        )

    # Calculate LCoE
    lcoe = investment / generation
    lcoe = lcoe.rename("lcoe")

    # Convert capacities and LCoE into an easier to handle format
    df = xr.merge(
        [
            capacities,
            lcoe,
            generation,
        ]
    ).to_dataframe()

    df = df.dropna().sort_values("lcoe")

    # Calculate annual generation across all sources + source classes
    df["cumulative generation"] = df.cumsum()["annual generation"]

    # Remove inf values (espc. for LCoE where generation = 0)
    df = df.replace([np.inf, -np.inf], np.nan).dropna()

    ## Reserve (fully and partially) capacities for domestic demand
    logger.info(f"Removing {demand:.0f} MWh of annual demand from potential supply.")

    # Helper column
    df["reserved share"] = np.nan

    # Supply class/technology which are fully consumed by domestic demand
    df.loc[df["cumulative generation"] <= demand, "reserved share"] = 1

    # Supply class/technology which can satisfy the remaining demand without being fully consumed
    idx = df.loc[df["reserved share"].isna()].index[0]
    _ds = df.loc[idx]

    # Before this technology/class
    residual_demand = demand - (_ds["cumulative generation"] - _ds["annual generation"])

    # Calculate the share of this technology reserved for domestic demand
    _ds["reserved share"] = (_ds["annual generation"] - residual_demand) / _ds[
        "annual generation"
    ]

    # Split class/technology into two entries:
    # one fully reserved and one with the remaining capacity
    _ds_remaining = _ds.copy()
    _ds_remaining["reserved share"] = 0

    _ds_reserved = _ds.copy()
    _ds_reserved["reserved share"] = 1

    _ds_remaining.loc[["capacities", "annual generation"]] *= 1 - _ds["reserved share"]
    _ds_reserved.loc[["capacities", "annual generation"]] *= _ds["reserved share"]

    _ds_reserved["cumulative generation"] -= _ds_remaining["annual generation"]

    # Remove the original class/technology entry and add the splitted ones
    df = df.drop(idx)
    df = pd.concat(
        [
            df,
            _ds_reserved.to_frame().transpose(),
            _ds_remaining.to_frame().transpose(),
        ]
    ).sort_values("lcoe")

    # Restore names for multiindex
    df.index.names = ["class", "technology"]

    # Indicate all remaining technologies/classes to be available for the ESC
    df["reserved share"] = df["reserved share"].fillna(0)

    # Save LCoEs as intermediary result
    df.to_csv(snakemake.output["lcoes"])

    # Exclude sources reserved for domestic demand
    df = df[df["reserved share"] == 0]
    if df.empty:
        logger.error(
            f"Insufficient renewable capacities ({generation.sum().item():.0f} MWh)"
            f" to satisfy domestic demand ({demand:.0f} MWh)."
        )

    df = df.reorder_levels(["technology", "class"]).sort_index()
    for idx, row in df.iterrows():
        # Select profiles as maximum dispatchable feed-in
        # Lower and upper clipping
        # max 100% (higher values due to EPS)
        # min 0.1% or 0. (easier for solving and non-relevant generation)
        p_max_pu = supply["profiles"].sel(
            {"technology": idx[0], "class": idx[1]}, drop=True
        )
        p_max_pu = p_max_pu.clip(min=0.0, max=1.0)
        p_max_pu = p_max_pu.where(lambda x: x > 0.0001, 0).to_numpy()

        if idx[0] in {"solar-utility", "offwind", "onwind"}:
            # These technologies directly feed-in electricity
            network.add(
                "Generator",
                name=f"{idx[0]} {idx[1]}",
                bus="electricity (exp)",
                p_nom_max=row["capacities"],
                p_nom_extendable=True,
                p_max_pu=p_max_pu,
                # Cost data for kW, here MW implicitly assumed
                capital_cost=1.0e3
                * calculate_annual_investment(idx[0], wacc, snakemake.input["costs"]),
                carrier=idx[0],
            )
        elif idx[0] == "csp-tower":
            # CSP is modelled with the generator as heat provider ("csp-tower #")
            # and additional components for heat storage ("csp-tower TES #"),
            # heat-to-electricity generator ("csp-tower power block #")
            # and a dedicated bus for each CSP plant class ("csp-tower bus #")
            bus_name = f"{idx[0]} {idx[1]} bus"
            network.add(
                "Bus",
                name=bus_name,
                carrier="heat",
                unit="MW",
            )

            name_ = f"{idx[0]}"
            network.add(
                "Generator",
                name=f"{name_} {idx[1]}",
                bus=bus_name,
                p_nom_max=row["capacities"],
                p_nom_extendable=True,
                p_max_pu=p_max_pu,
                # Cost data for kW, here MW implicitly assumed
                capital_cost=1.0e3
                * calculate_annual_investment(name_, wacc, snakemake.input["costs"]),
                carrier=f"{idx[0]}",
            )

            name_ = f"{idx[0]} TES"
            network.add(
                "Store",
                name=f"{name_} {idx[1]}",
                bus=bus_name,
                e_nom_extendable=True,
                e_cyclic=True,
                # Scale from kW_e in cost data to MW_e
                capital_cost=1e3
                * calculate_annual_investment(name_, wacc, snakemake.input["costs"]),
            )

            name_ = f"{idx[0]} power block"
            efficiency_ = (
                read_efficiencies(
                    snakemake.input["efficiencies"], snakemake.wildcards["year"]
                )
                .query("process == @name_")["efficiency"]
                .item()
            )
            network.add(
                "Link",
                name=f"{name_} {idx[1]}",
                bus0=bus_name,
                bus1="electricity (exp)",
                p_nom_extendable=True,
                efficiency=efficiency_,
                # capital_cost are scaled with efficiency
                # Scale from kW_e in cost data to MW_e
                capital_cost=1e3
                * efficiency_
                * calculate_annual_investment(name_, wacc, snakemake.input["costs"]),
            )

        else:
            logger.error(f"Unknown technology '{idx[0]}'.")

    network.export_to_netcdf(snakemake.output["network"])

    supply.close()

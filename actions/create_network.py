# SPDX-FileCopyrightText: 2020-2023 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later


import logging
import pickle
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pypsa
from _helpers import (
    calculate_annual_investment,
    calculate_annuity,
    configure_logging,
    extract_technology,
    get_bus_unit,
    read_efficiencies,
)

logger = logging.getLogger(__name__)


def _do_units_match(unit1, unit2):
    """Rough checker if two units match.

    Catches mismatches between e.g. t, m^3 and W.
    Does not catch order of magnitude mismatches by SI prefixes."""

    def _stripper(u):
        import re

        # remove optional Si prefixes 'M', 'k' and per hour ('/h', 'h') indicators
        u = re.match("[Mk]?(.+?)\/?h?$", re.split("_|-", u)[0]).groups()[0]

        # Specific indicators to remove, a bit hacky and may cause problems they do not match. Hotfix.
        for p in ["CO2", "/km"]:
            u = u.replace(p, "")

        return u

    return _stripper(unit1) == _stripper(unit2)


def create_network():
    """Create the pypsa network scaffolding for the ESC."""

    # Modify PyPSA 'Link' component to allow for 2 output busses by overwriting component_attrs
    # c.f. https://www.pypsa.org/examples/chp-fixed-heat-power-ratio.html

    # Load additional components
    with open(snakemake.input["additional_components"], "rb") as f:
        override_component_attrs = pickle.load(f)

    # Create network with modified link-component
    network = pypsa.Network(override_component_attrs=override_component_attrs)

    # Load network components from csv files
    network.import_from_csv_folder(snakemake.input["network"])

    # Equally weighted snapshots, year defined via config
    year = int(snakemake.params["era_year"])
    snapshots = pd.date_range(str(year), str(year + 1), freq="H", inclusive="left")
    network.set_snapshots(snapshots)

    return network


def attach_efficiencies(network):
    """Attach dedicated efficiencies from file.

    The efficiencies are from an additional csv file and added to the links in the pypsa network
    Format for efficiencies.csv file:
    * "from" and "to" must substrings of the bus names
    * "process" must be a substring of the name of the link

    Return
    ------
    network : pypsa.network
        network with external efficiencies attached to all links.

    """
    efficiencies = read_efficiencies(
        snakemake.input["efficiencies"], snakemake.wildcards["year"]
    )

    def get_efficiency(tech, src_bus, tar_bus):
        tech = extract_technology(tech)
        src = extract_technology(src_bus)
        tar = extract_technology(tar_bus)

        efficiency = efficiencies[
            (efficiencies["process"] == tech)
            & (efficiencies["to"] == tar)
            & (efficiencies["from"] == src)
        ]

        if efficiency.empty is True:
            return np.nan

        # Check if all units match
        src_unit = get_bus_unit(src_bus, network)
        tar_unit = get_bus_unit(tar_bus, network)
        unit_mismatch = None
        if (efficiency[["from_unit", "to_unit"]] == "p.u.").any(
            axis=None
        ) == True:  #'==' b/c pd.any returns np.bool
            if (efficiency[["from_unit", "to_unit"]] == "p.u.").all(
                axis=None
            ) == False:  #'==' b/c pd.any returns np.bool
                unit_mismatch = "One efficiency in [p.u.], but the other one not."
            elif _do_units_match(src_unit, tar_unit) is False:
                unit_mismatch = f"Unit of bus {src_bus} [{src_unit}] does not match {tar_bus} [{tar_unit}]."

        elif _do_units_match(src_unit, efficiency["from_unit"].item()) is False:
            unit_mismatch = (
                f"Source bus unit {src_bus} [{src_unit}] does not match unit "
                f'in registered efficiencies.csv [{efficiency["from_unit"].item()}].'
            )
        elif _do_units_match(tar_unit, efficiency["to_unit"].item()) is False:
            unit_mismatch = (
                f"Target bus unit {tar_bus} [{tar_unit}] does not match unit "
                f'in registered efficiencies.csv [{efficiency["to_unit"].item()}].'
            )

        if unit_mismatch:
            raise ValueError(f"Mismatching units for {tech}: {unit_mismatch}.")

        return efficiency["efficiency"].item()

    links = network.links

    for idx, row in links.iterrows():
        lead_efficiency = get_efficiency(
            extract_technology(row.name), row["bus0"], row["bus1"]
        )
        if np.isnan(lead_efficiency):
            logger.error(
                f"No efficiency found for link {row.name} between "
                f"{row['bus0']} and {row['bus1']}."
            )
        else:
            links.loc[idx, "efficiency"] = lead_efficiency

        additional_buses = {
            c for c in links.columns if c.startswith("bus") and row[c] != ""
        } - {"bus0", "bus1"}
        for b in additional_buses:
            # by design decision all buses busn (n>1, e.g. bus2, bus3, ...) either:
            # case 1. contribute to the output to bus1, e.g. bus2 feeds into bus1
            # or
            # case 2. are few by bus0.
            # Try to retrieve both efficiencies: One should always be nan, the other one is taken.
            # If either one are nan or not nan, throw an error.
            #
            # Case 1.:
            # Efficiencies are provided for the conversion from bus2 to bus1
            # and are thus weighted by the primary efficiency of bus1
            # Efficiencies are have to become negative to correctly account for the flow.
            follow_efficiency = (
                (-1)
                * lead_efficiency
                / get_efficiency(extract_technology(row.name), row[b], row["bus1"])
            )

            # Case 2.:
            # Efficiencies are provided for conversion from bus0 to busn.
            # This type of efficiency does not need to be adjusted.
            regular_efficiency = efficiency = get_efficiency(
                extract_technology(row.name), row["bus0"], row[b]
            )

            e = np.array([regular_efficiency, follow_efficiency])

            if np.isnan(e).all():
                logger.error(
                    f"No efficiency found for link {row.name} between "
                    f"{row[b]} <-> ({row['bus0']} or {row['bus1']})."
                )
            elif np.isnan(e).sum() == 0:
                logger.error(
                    f"Two efficiencies found for link {row.name} between "
                    f"{row[b]} <-> ({row['bus0']} and {row['bus1']}). "
                    f"Efficiency must by unambigious."
                )
            else:
                links.loc[idx, "efficiency" + b.replace("bus", "")] = e[~np.isnan(e)][0]

    return network


def override_costs_for_special_cases(n):
    # battery inverter represented by two links (charging and discharging),
    # while costs in cost data are for bidirectional inverter --> correction here
    links = n.links
    idx = links.filter(like="battery inverter", axis=0).index
    links.loc[idx, "capital_cost"] /= 2.0

    return n


def attach_costs(network):
    """
    Attach the overnight investment costs (capital costs) to the network.

    Costs are calculated from investment costs and FOM using EAC method
    and wacc as specified via config/snakemake.input files.
    Components name need to follow the scheme '<name> (<exp|imp>)'
    where '<name>' must correspond to the component in the costs.csv file.

    Requires efficiencies of the network to be already attached.
    """
    wacc = pd.read_csv(
        snakemake.input["wacc"], comment="#", index_col="region", keep_default_na=False
    )
    wacc = wacc.loc[snakemake.wildcards["from"], scenario["wacc"]]
    wacc *= scenario["modifiers"]["wacc"]  # Apply scenario modifier

    costs = pd.read_csv(snakemake.input["costs"], index_col=["technology", "parameter"])

    def attach_component_costs(network, component):
        components = getattr(network, component)
        for idx, row in components.iterrows():
            try:
                tech = costs.loc[extract_technology(idx)]
            except KeyError:
                logger.warning(f"No cost assumptions found for {idx}.")
                continue

            if "commodity" in tech.index:
                logger.info(f"Adding commodity costs for {tech.index}")
                assert _do_units_match(
                    get_bus_unit(row["bus"], network),
                    tech.loc["commodity"]["unit"].replace(
                        "EUR/", ""
                    ),  # Assume unit for fuel/commodities is EUR/<bus unit>
                ), "Mismatching units"
                components.loc[idx, "marginal_cost"] = tech.loc["commodity"]["value"]
                continue

            # Compare units between bus and cost data; scale investment on-demand
            investment_unit = (
                tech.loc["investment"]["unit"]
                .replace("EUR/", "")
                .replace("(", "")
                .replace(")", "")
            )

            investment_factor = 1.0

            # Determine unit of the bus - depends on component type
            # stores -> attached to one bus, situation clear
            # generators -> attached to one bus, situation clear
            # links -> attached to >=2 buses, use additional column in .csv to determine the bus to scale to
            bus_unit = bus0_unit = bus1_unit = None
            if component in ["stores", "generators"]:
                bus_unit = network.buses.loc[row["bus"]]["unit"]
            elif component == "links":
                if row["scale_costs_based_on"] not in ["bus0", "bus1"]:
                    raise NotImplementedError(
                        f"Scaling of costs for link {idx} to others than bus0 or bus1 not implemented."
                    )

                bus_unit = network.buses.loc[row[row["scale_costs_based_on"]]]["unit"]
                if row["scale_costs_based_on"] == "bus1":
                    # bus1 is output by convention; cost in output units, i.e. scale to input units
                    investment_factor *= row["efficiency"]

            # Consistency check: Correct units (ignore first letter for prefixed values)
            # Warning: Does not catch case were bus unit is only a single letter
            if _do_units_match(investment_unit, bus_unit) is False:
                raise ValueError(
                    f'Could not find matching cost data for {component} "{idx}": '
                    f"Expected {bus_unit} based on network, but found {investment_unit} in cost data."
                )

            prefix_bus_unit = bus_unit[0]
            prefix_investment_unit = investment_unit[0]

            if prefix_bus_unit == prefix_investment_unit:
                investment_factor *= 1.0
            elif prefix_bus_unit == "M" and prefix_investment_unit == "k":
                investment_factor *= 1.0e3
            else:
                raise ValueError(
                    f"Cannot scale between {prefix_bus_unit} and {prefix_investment_unit} "
                    f"for {idx} costs."
                )

            # Some technologies are without FOM values
            # (e.g. battery capacity where FOM is attributed to the link/inverter/charger capacities)
            try:
                fom = tech.loc["FOM", "value"]
            except KeyError:
                logger.info(f"No FOM for {idx}, assuming 0%.")
                fom = 0.0

            logger.info(
                f"Rescaling investment for {idx} by {investment_factor} to match units and efficiencies."
            )
            capital_cost = calculate_annuity(
                tech.loc["investment", "value"] * investment_factor,
                fom,
                tech.loc["lifetime", "value"],
                wacc,
            )
            components.loc[idx, "capital_cost"] = capital_cost

        return network

    network = attach_component_costs(network, "links")
    network = attach_component_costs(network, "stores")
    network = attach_component_costs(network, "generators")

    return network


def scale_transportation_with_distance(
    n, link_types=["HVDC overhead", "HVDC submarine", "pipeline", "submarine pipeline"]
):
    """Scales the cost and efficiency of specific links (transport options) by distance.

    Only implemented for:
     * HVDC overhead|submarine
     * H2 (g) [submarine] pipeline
     * CH4 (g) [submarine] pipeline
    Names of the components must match the regex, otherwise no scaling happens.
    Check the 'length' attribute of a link to see if it has been scaled and to which distance.
    A distance of 0 is forcefully set to distance=1 for compressors to keep this component
    and conversion in the supply chain including its connection to the electricity bus.
    Costs and efficiencies are - if affected - scaled accordingly.
    Submarine distances are substracted from HVDC overhead and onland pipelines from their distance.
    """

    links = n.links
    distances = pd.read_csv(
        snakemake.input["distances"], comment="#", quotechar='"', keep_default_na=False
    )

    # These are the components to scale
    # they may behave differently, case-by-case decision below
    mapping = {
        "HVDC overhead": {
            "name_pattern": "HVDC overhead$",  # pattern to select components from n.links using regex
            "distance_type": "as-the-crow-flies",  # entry as in snakemake.input['distances']
            "reduce_by_distance_type": "as-the-hake-swims",  # entry as in snakemake.input['distances']
            "detour_factor_key": "transmission_line",  # from snakemake.config
        },
        "HVDC submarine": {
            "name_pattern": "HVDC submarine$",
            "distance_type": "as-the-hake-swims",
            "detour_factor_key": "transmission_line",
        },
        "pipeline": {
            "name_pattern": "(H2|CH4) \(g\) pipeline$",
            "distance_type": "as-the-crow-flies",
            "reduce_by_distance_type": "as-the-hake-swims",
            "detour_factor_key": "pipeline",
        },
        "submarine pipeline": {
            "name_pattern": "(H2|CH4) \(g\) submarine pipeline$",
            "distance_type": "as-the-hake-swims",
            "detour_factor_key": "pipeline",
        },
        "pipeline compressor": {
            "name_pattern": "^(H2|CH4) \(g\) pipeline compressor(\s\(exp|imp\))?",
            "distance_type": "as-the-crow-flies",
            "detour_factor_key": "pipeline",
        },
    }

    distances = distances[
        (distances["region_a"] == snakemake.wildcards["from"])
        & (distances["region_b"] == snakemake.wildcards["to"])
    ].set_index("type")["value"]

    efs = [c for c in links.columns if c.startswith("efficiency")]

    for link_type in link_types:
        # Relevant properties for this link_type
        m = mapping[link_type]

        # All links in the network matching the link_type
        idx = links.filter(regex=m["name_pattern"], axis=0).index
        if idx.empty is True:
            continue  # Nothing to do

        distance = distances[m["distance_type"]]
        if m.get("reduce_by_distance_type", False):
            distance -= distances[m["reduce_by_distance_type"]]
        distance *= snakemake.config["detour_factors"][m["detour_factor_key"]]

        # For bookkeeping
        links.loc[idx, "length"] = distance

        # capital_cost in EUR/km, scale linearly
        links.loc[idx, "capital_cost"] *= distance

        # efficiencies which scale by power law per 1000km
        links.loc[idx, "efficiency"] = (links.loc[idx, "efficiency"]) ** (
            distance / 1.0e3
        )

        if links.loc[idx, efs].notnull().to_numpy().sum() != 1:
            raise NotImplementedError(
                f"Scaling of link {idx} failed. "
                f"Expected exactly one efficiency, but found more or less."
            )

    return n


def add_shipping(n):
    """Adds optional shipping routes to the network.

    Checks whether file "<ESC>/ships.csv" exists and - if it does - constructs
    a shipping route with multiple convoys (as optimisation options) for this route
    using standard PyPSA components.

    Limitation: ONLY SUPPORTS ONE SHIPPING CONNECTION PER NETWORK AT THE MOMENT!
    """

    fn = Path(snakemake.input["network"]) / "ships.csv"

    # Network without shipping routes
    if not fn.exists():
        return n
    else:
        ships = pd.read_csv(fn, comment="#", index_col="name")

    props = pd.read_csv(
        snakemake.input["shipping_properties"],
        comment="#",
        index_col=["name", "variable"],
    )

    distances = pd.read_csv(
        snakemake.input["distances"], comment="#", quotechar='"', keep_default_na=False
    )

    wacc = pd.read_csv(
        snakemake.input["wacc"], comment="#", index_col="region", keep_default_na=False
    )
    wacc = wacc.loc[snakemake.wildcards["from"], scenario["wacc"]]
    wacc *= scenario["modifiers"]["wacc"]  # Apply scenario modifier

    costs = pd.read_csv(snakemake.input["costs"], index_col=["technology", "parameter"])

    if len(ships.index) != 1 and n.name != "LOHC shipping":
        logger.warning(
            "More than one entry in ships.csv only supported for LOHC shipping."
        )

    ship = ships.iloc[0]

    # Which fuel to use for ship propulsion:
    # Either 'cargo' (default) or name of a shipping fuel in the cost csv file
    assert (
        "use_fuel" in ship
    ), "use_fuel needs to be defined in ships.csv (either shipping fuel or 'cargo')"
    fuel_to_use = ship["use_fuel"]

    props = props.loc[ship.name]

    loading_time = int(np.floor(props.loc["(un-) loading time", "value"]))
    unloading_time = loading_time
    loading_rate_pu = 1.0 / loading_time
    unloading_rate_pu = loading_rate_pu

    distance = distances.query(
        f"""region_a == '{snakemake.wildcards['from']}' and """
        f"""region_b == '{snakemake.wildcards['to']}' and """
        f"""type == 'sea route'""",
        engine="python",
    )["value"].item()

    travel_time = int(np.ceil(distance / props.loc["average speed", "value"]))

    # Round trip time for a convoy (loading, travel, unloading, return trip)
    round_trip_time = loading_time + travel_time + unloading_time + travel_time

    # Number of full journeys (round-trip-journey) possible for convoy along sea route
    journeys = int(np.floor(n.snapshots.shape[0] / round_trip_time))

    # By constructing the tightest shipping schedule starting at the beginning of the year
    # we have this amount of hours were the importing habour is not served...
    annual_shipping_gap = n.snapshots.shape[0] % round_trip_time

    # ... as we later construct additional shipping convoys by simply shifting the schedule,
    # this will create a nasty gap in the supply chain, resulting in weired results in the optimisation.
    # We avoid this by smoothing the supply: the shipping duration is artificially prolonged to reduce this gap
    # Can be thought of as something like a buffer, which is near identically distributed across all journeys
    additional_forward_travel_time = int(np.floor(annual_shipping_gap / journeys / 2))
    # return trip can take a bit longer (max 1 additional snapshot)
    additional_return_travel_time = int(
        np.floor(annual_shipping_gap / journeys - additional_forward_travel_time)
    )

    forward_travel_time = travel_time + additional_forward_travel_time
    return_travel_time = travel_time + additional_return_travel_time

    updated_round_trip_time = (
        loading_time + forward_travel_time + unloading_time + return_travel_time
    )
    logger.info(
        f"Increasing the round-trip travel time from {round_trip_time}h to "
        f"{updated_round_trip_time}h (+{(updated_round_trip_time/round_trip_time-1)*100:.2f}%) "
        f"to achieve more levelled supply by ship."
    )

    # One round-trip loading schedule for earliest convoy in year
    loading_schedule = np.concatenate(
        (
            [loading_rate_pu] * loading_time,
            [0] * forward_travel_time,
            [0] * unloading_time,
            [0] * return_travel_time,
        )
    )

    # One round-trip unloading schedule for earliest convoy in year
    unloading_schedule = np.concatenate(
        (
            [0] * loading_time,
            [0] * forward_travel_time,
            [unloading_rate_pu] * unloading_time,
            [0] * return_travel_time,
        )
    )

    # Numbers of convoys (base convoy + convoys which can be loaded without competing for the
    #  loading infrastructure while the base convoy is on its way)
    # if loading_time == unloading_time there this approach results in no clashes for the unloading infrastruct.
    convoy_number = 1 + int(
        np.floor(
            (forward_travel_time + return_travel_time + unloading_time) / loading_time
        )
    )

    shift = (forward_travel_time + unloading_time + return_travel_time) % unloading_time
    shift_per_convoy = int(np.floor(shift / convoy_number))

    # Create full year schedule for loading:
    # Left-over days at end of year (which can not be used for a full round-trip journey)
    # are filled with 0s (=no journey/anchored)
    loading_schedule = np.concatenate([loading_schedule] * journeys)
    tmp = np.zeros(n.snapshots.shape[0])
    tmp[: loading_schedule.shape[0]] = loading_schedule
    loading_schedule = tmp

    # Create full year schedule for unloading:
    # (Basically the same as for loading, could use np.roll here as rates and durations for loading
    #  and unloading are identical in the current model version)
    # Left-over days at end of year (which can not be used for a full round-trip journey)
    # are filled with 0s (=no journey/anchored)
    unloading_schedule = np.concatenate([unloading_schedule] * journeys)
    tmp = np.zeros(n.snapshots.shape[0])
    tmp[: unloading_schedule.shape[0]] = unloading_schedule
    unloading_schedule = tmp

    ## How the schedules look like
    # plt.plot(loading_schedule, label='loading')
    # plt.plot(unloading_schedule, label='unloading')
    # plt.legend()

    # Calculate energy transport efficiency for the trip

    # Boil-off losses are only considered for the outward journey.
    # technically the cargo hold should be kept at cryo temperatures and uncontaminated also during
    # the inward journey. Neglect boil-off during inward journey here, as this can be expected to be
    # significantly smaller than the onboard energy demand (all boil-off is certainly consumed by energy demand)
    boil_off = (1 - props.loc["boil-off", "value"] / 100) ** forward_travel_time

    # Energy demand in MWh for ship propulsion (outward and return journey)
    energy_demand = (
        2
        * distance
        * props.loc["energy demand", "value"]
        / props.loc["capacity", "value"]
    )

    if fuel_to_use == "cargo":
        # take whatever requires more energy (boil-off can be used by propulsion or propulsion uses cargo)
        shipping_efficiency = np.min([1 - energy_demand, boil_off])
        assert shipping_efficiency > 0, "Shipping (Link) efficiency must be > 0"
        logging.info("Using shipping cargo as fuel source.")
    else:
        shipping_efficiency = boil_off

        logger.info(
            f"""
            Using external shipping fuel: {fuel_to_use}
            Adding necessary network components (fuel bus and generator for purchasing fuel).
            """
        )
        n.add(
            "Bus",
            name=f"{fuel_to_use} bus",
            carrier=f"{fuel_to_use}",
            unit="MWh",
        )

        n.add(
            "Generator",
            name=f"{fuel_to_use} purchase",
            bus=f"{fuel_to_use} bus",
            carrier=f"{fuel_to_use}",
            capital_cost=0.0,
            marginal_cost=costs.loc[fuel_to_use, "fuel"]["value"],
            p_nom_extendable=True,
            unit="MWh",
        )

    # Additional energy losses from (un-) loading the cargo
    loading_efficiency = 1 - props.loc["(un-) loading losses", "value"] / 100
    unloading_efficiency = loading_efficiency

    # Calculate investment costs per gross MWh capacity
    costs = costs.loc[ship.name]

    # Consistency check: whether units match
    unit_costs = costs.loc["capacity"]["unit"]
    unit_bus = network.buses.loc[ship["bus0"]]["unit"]
    if _do_units_match(unit_bus, unit_costs) is False:
        raise ValueError(
            f"Unit mismatch for {ship.name} between network ({unit_bus}) "
            f"and cost database ({unit_costs})."
        )

    try:
        capital_cost = calculate_annuity(
            costs.loc["investment", "value"],
            costs.loc["FOM", "value"],
            costs.loc["lifetime", "value"],
            wacc,
        )
        capital_cost /= costs.loc["capacity", "value"]
    except:
        raise ValueError(
            f"Exception calculating capital cost for {ship.name}."
            f"Missing cost or shipping property entries"
        )

    if distance == 0:
        # Treat the special case, where distance between exporter and importer is zero.
        # (e.g. same exporter as importer region)

        logger.info(
            'No distance between exporter and importer. Adding a direct pseudo connection without shipping schedule.'
        )

        convoy_number = 1

        loading_schedule = np.ones(n.snapshots.shape[0])
        unloading_schedule = loading_schedule

    else:
        logger.info(f"Adding {convoy_number} shipping convoys to shipping route.")

    for i in range(convoy_number):
        ship_bus = f"{ship.name} convoy {i+1}"

        n.add(
            "Bus",
            name=f"{ship_bus} (exp)",
            carrier=network.buses.loc[ship.loc["bus0"], "carrier"],
            unit=network.buses.loc[ship.loc["bus0"], "unit"],
        )

        n.add(
            "Link",
            name=f"{ship_bus} loading",
            bus0=ship.loc["bus0"],
            bus1=f"{ship_bus} (exp)",
            efficiency=loading_efficiency,
            # Capacity expansion at point of export/depature
            # ship capacity taken into account as gross capacity before transport losses/demand
            p_nom_extendable=True,
            # Loading extracts energy from bus and happens at max rate and at fixed times
            # Rolling the schedules ensures there is no overlap between convoys
            # p_min_pu=np.roll(
            #    loading_schedule, i * (loading_time + shift_per_convoy)
            # ),  # np.zeros_like(loading_schedule),
            p_min_pu=0.0,
            p_max_pu=np.roll(loading_schedule, i * (loading_time + shift_per_convoy)),
        )

        n.add(
            "Store",
            name=f"{ship_bus} cargo (exp)",
            bus=f"{ship_bus} (exp)",
            e_nom_extendable=True,
            # Ships may starting at end of year and start deliver at the beginning of the year
            e_cyclic=True,
            # Full capital cost of the ship
            capital_cost=capital_cost,
        )

        n.add(
            "Bus",
            name=f"{ship_bus} (imp)",
            carrier=network.buses.loc[ship.loc["bus0"], "carrier"],
            unit=network.buses.loc[ship.loc["bus0"], "unit"],
        )

        n.add(
            "Link",
            name=f"{ship_bus} unloading",
            bus0=f"{ship_bus} (imp)",
            bus1=ship.loc["bus1"],
            efficiency=unloading_efficiency,
            p_nom_extendable=True,
            # Unloading at max rate and at fixed times
            # Rolling the schedules ensures there is no overlap between convoys
            # p_min_pu=np.roll(
            #    unloading_schedule, i * (unloading_time + shift_per_convoy)
            # ),  # np.zeros_like(unloading_schedule),
            p_min_pu=0.0,
            p_max_pu=np.roll(
                unloading_schedule, i * (unloading_time + shift_per_convoy)
            ),
        )

        n.add(
            "Link",
            name=f"{ship_bus} trip demand & losses",
            bus0=f"{ship_bus} (exp)",
            bus1=f"{ship_bus} (imp)",
            efficiency=shipping_efficiency,
            p_nom_extendable=True,
            # If not cargo but external fuel is used: Connect to fuel supplying bus
            bus2=f"{fuel_to_use} bus" if fuel_to_use != "cargo" else "",
            efficiency2=(-1) * energy_demand if fuel_to_use != "cargo" else np.nan,
            # efficiency2=(-1)*energy_demand/shipping_efficiency if fuel_to_use != "cargo" else "",
        )

        # Special case LOHC:
        # requires a return channel for unloaded LOHC using the same ship and thus a dedicated store component
        # Warning: Partially hardcoded!
        # Note: Relies on information from ships.csv
        # Note: Linked to a extra_functionality in solve_network.py.ipynb fixing the two store's capacities
        if n.name == "LOHC shipping":
            if fuel_to_use != "cargo":
                raise NotImplementedError(
                    "Use of external fuel (instead of cargo) is not yet implemented for LOHC shipping. Efficiency of LOHC consuming ship-internal links need to be adjusted."
                )

            n.add(
                "Bus",
                name=f"{ship_bus} LOHC (used)",
                carrier=n.buses.loc["LOHC (used) (exp)", "carrier"],
                unit=n.buses.loc["LOHC (used) (exp)", "unit"],
            )

            n.add(
                "Store",
                name=f"{ship_bus} cargo LOHC (used)",
                bus=f"{ship_bus} LOHC (used)",
                e_nom_extendable=True,
                e_cyclic=True,
                # Capital costs already considered by the cargo store for loaded LOHC
                capital_cost=0,
            )

            n.add(
                "Link",
                name=f"{ship_bus} LOHC (used) loading",
                bus0=ships.loc["LOHC (used) transport ship", "bus0"],
                bus1=f"{ship_bus} LOHC (used)",
                efficiency=1.0,
                p_nom_extendable=True,
                # Loading extracts energy from bus and happens at max rate and at fixed times
                # Rolling the schedules ensures there is no overlap between convoys
                p_min_pu=np.roll(
                    unloading_schedule, i * (unloading_time + shift_per_convoy)
                ),  # np.zeros_like(unloading_schedule),
                p_max_pu=np.roll(
                    unloading_schedule, i * (unloading_time + shift_per_convoy)
                ),
            )

            n.add(
                "Link",
                name=f"{ship_bus} LOHC (used) unloading",
                bus0=f"{ship_bus} LOHC (used)",
                bus1=ships.loc["LOHC (used) transport ship", "bus1"],
                efficiency=1.0,
                p_nom_extendable=True,
                # Unloading at max rate and at fixed times
                # Rolling the schedules ensures there is no overlap between convoys
                p_min_pu=np.roll(
                    loading_schedule, i * (loading_time + shift_per_convoy)
                ),  # np.zeros_like(loading_schedule),
                p_max_pu=np.roll(
                    loading_schedule, i * (loading_time + shift_per_convoy)
                ),
            )

            # Propulsion energy demand
            # if LOHc cargo is to be used for propulsion, make the link a multi-link
            # such that unloaded LOHC is accounted for
            if fuel_to_use == "cargo":
                # Ratio LOHC / H2 in loaded LOHC
                loaded_unloaded_ratio = read_efficiencies(
                    snakemake.input["efficiencies"], snakemake.wildcards["year"]
                )
                loaded_unloaded_ratio = loaded_unloaded_ratio.loc[
                    (loaded_unloaded_ratio["process"] == "LOHC dehydrogenation")
                    & (loaded_unloaded_ratio["from"] == "LOHC (loaded)")
                    & (loaded_unloaded_ratio["to"] == "LOHC (used)")
                ]["efficiency"].item()
                n.links.loc[
                    f"{ship_bus} trip demand & losses", "bus2"
                ] = f"{ship_bus} LOHC (used)"
                n.links.loc[f"{ship_bus} trip demand & losses", "efficiency2"] = (
                    1 - shipping_efficiency
                ) * loaded_unloaded_ratio
    return n


# In[ ]:


if __name__ == "__main__":
    scenario = snakemake.params["scenario"]

    configure_logging(snakemake)

    # DO NOT CHANGE THIS ORDER; dependencies between methods not explicit
    network = create_network()
    network = attach_efficiencies(network)
    network = attach_costs(network)

    network = scale_transportation_with_distance(network)

    network = add_shipping(network)

    network = override_costs_for_special_cases(network)

    network.export_to_netcdf(snakemake.output["network"])

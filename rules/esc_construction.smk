# SPDX-FileCopyrightText: 2020-2021 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later


rule create_additional_components:
    output:
        additional_components="resources/additional_components.pkl",
    threads: 1
    log:
        python="logs/create_additional_components.log",
    notebook:
        "../actions/create_additional_components.py.ipynb"


rule create_network:
    input:
        efficiencies="data/efficiencies.csv",
        costs=technology_data("../technology-data/outputs/costs_{year}.csv"),
        wacc="data/wacc.csv",
        distances="data/distances.csv",
        shipping_properties="data/shipping.csv",
        network="escs/{esc}",
        additional_components="resources/additional_components.pkl",
    output:
        network="resources/networks/{scenario}/{year}/{esc}/{from}-{to}/network.nc",
    threads: 1
    log:
        python="logs/{scenario}/{year}/{esc}/{from}-{to}/create_network.log",
        notebook="logs/{scenario}/{year}/{esc}/{from}-{to}/create_network.ipynb",
    notebook:
        "../actions/create_network.py.ipynb"


def demand_file(wildcards):
    # Allow for custom overwrite of annual electricity demand for exporters
    demand_d = {
        "gegis": "resources/demand_TRACES_2013.nc",
        "custom": "data/overwrite/demand.csv",
    }

    choice = config["scenarios"][wildcards.scenario].get(
        "synthetic_demand", config["scenarios"]["default"]["synthetic_demand"]
    )
    return demand_d[choice.lower()]


rule attach_supply:
    input:
        supply="resources/supply_TRACES_2013.nc",
        demand=demand_file,
        costs=technology_data("../technology-data/outputs/costs_{year}.csv"),
        wacc="data/wacc.csv",
        network="resources/networks/{scenario}/{year}/{esc}/{from}-{to}/network.nc",
        additional_components="resources/additional_components.pkl",
    output:
        network=(
            "resources/networks_supplied/{scenario}/{year}/{esc}/{from}-{to}/network.nc"
        ),
        lcoes=(
            "resources/networks_supplied/{scenario}/{year}/{esc}/{from}-{to}/lcoes.csv"
        ),
    threads: 1
    log:
        python="logs/{scenario}/{year}/{esc}/{from}-{to}/attach_supply.log",
        notebook="logs/{scenario}/{year}/{esc}/{from}-{to}/attach_supply.ipynb",
    notebook:
        "../actions/attach_supply.py.ipynb"

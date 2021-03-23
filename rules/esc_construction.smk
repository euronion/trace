# SPDX-FileCopyrightText: 2020-2021 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later

rule create_additional_components:
    output:
        additional_components="resources/additional_components.pkl"
    threads: 1
    log:
        python="logs/create_additional_components.log",
    notebook:
        "../actions/create_additional_components.py.ipynb"

rule create_network:
    input:
        efficiencies="data/efficiencies.csv",
        costs=f"../technology-data/outputs/costs_{config['scenario']['year']}.csv",
        wacc="data/wacc.csv",
        distances="data/distances.csv",
        shipping_properties="data/shipping.csv",
        network="escs/{esc}",
        additional_components="resources/additional_components.pkl"
    output:
        network="resources/networks/"+SCENARIO_FOLDER+"/{esc}/{from}-{to}/network.nc"
    threads: 1
    log:
        python="logs/"+SCENARIO_FOLDER+"/create_network/{esc}/{from}-{to}.log",
        notebook="logs/"+SCENARIO_FOLDER+"/create_network/{esc}/{from}-{to}.ipynb"
    notebook:
        "../actions/create_network.py.ipynb"


# Allow for custom overwrite of annual electricity demand for exporters
demand_d = {
    "gegis": "resources/demand_TRACES_2013.nc",
    "custom": "data/overwrite/demand.csv"
}
demand_i = demand_d[config["scenario"]["synthetic_demand"].lower()]

rule attach_supply:
    input:
        supply="resources/supply_TRACES_2013.nc",
        demand=demand_i,
        costs=f"../technology-data/outputs/costs_{config['scenario']['year']}.csv",
        wacc="data/wacc.csv",
        network="resources/networks/"+SCENARIO_FOLDER+"/{esc}/{from}-{to}/network.nc",
        additional_components="resources/additional_components.pkl"
    output:
        network="resources/networks_supplied/"+SCENARIO_FOLDER+"/{esc}/{from}-{to}/network.nc",
        lcoes="resources/networks_supplied/"+SCENARIO_FOLDER+"/{esc}/{from}-{to}/lcoes.csv"
    threads: 1
    log:
        python="logs/"+SCENARIO_FOLDER+"/attach_supply/{esc}/{from}-{to}.log",
        notebook="logs/"+SCENARIO_FOLDER+"/attach_supply/{esc}/{from}-{to}.ipynb"
    notebook:
        "../actions/attach_supply.py.ipynb"

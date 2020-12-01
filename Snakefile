configfile: "config.yaml"

rule create_network:
    input:
        efficiencies="data/efficiencies.csv",
        network="escs/{esc}"
    output:
        directory("resources/networks/{esc}/{from}")
    log:
        notebook="logs/notebooks_processed/create_network/{esc}/{from}.ipynb"
    notebook:
        "actions/create_network.py.ipynb"


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
        costs="data/costs.csv",
        wacc="data/wacc.csv",
        network="resources/networks/{esc}/{from}"
    output:
        network=directory("resources/networks_supplied/{esc}/{from}")
    log:
        notebook="logs/notebooks_processed/attach_supply/{esc}/{from}.ipynb"
    notebook:
        "actions/attach_supply.py.ipynb"
        
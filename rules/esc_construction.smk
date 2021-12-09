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
    params:
        scenario=lambda w: get_scenario(w["scenario"]),
        era_year=config["GlobalEnergyGIS"]["era_year"],
    log:
        python="logs/{scenario}/{year}/{esc}/{from}-{to}/create_network.log",
        notebook="logs/{scenario}/{year}/{esc}/{from}-{to}/create_network.ipynb",
    notebook:
        "../actions/create_network.py.ipynb"


def get_import_profile_path(wildcards):
    """Return the import_profile filepath for the wildcards scenario or the default profile."""
    ip = config["scenarios"][wildcards.scenario].get(
        "import_profile", config["scenarios"]["default"]["import_profile"]
    )

    return f"data/import_profiles/{ip}.csv"


rule attach_import_profile:
    message:
        "Attaching import profile ('ip') to network."
    input:
        network="resources/networks/{scenario}/{year}/{esc}/{from}-{to}/network.nc",
        additional_components="resources/additional_components.pkl",
        import_profile=get_import_profile_path,
    output:
        network="resources/networks_ip/{scenario}/{year}/{esc}/{from}-{to}/network.nc",
    threads: 1
    params:
        scenario=lambda w: get_scenario(w["scenario"]),
    log:
        python="logs/{scenario}/{year}/{esc}/{from}-{to}/attach_import_profile.log",
        notebook=(
            "logs/{scenario}/{year}/{esc}/{from}-{to}/attach_import_profile.ipynb"
        ),
    notebook:
        "../actions/attach_import_profile.py.ipynb"


def demand_file(wildcards):
    # Allow for custom overwrite of annual electricity demand for exporters
    demand_d = {
        "gegis": "resources/demand_annual_TRACES_2013.csv",
        "custom": "data/overwrite/demand.csv",
    }

    choice = config["scenarios"][wildcards.scenario].get(
        "synthetic_demand", config["scenarios"]["default"]["synthetic_demand"]
    )
    return demand_d[choice.lower()]


rule attach_supply:
    message:
        "Attaching RES supply ('as') to network."
    input:
        supply="resources/supply_TRACES_{era_year}.nc".format(
            era_year=config["GlobalEnergyGIS"]["era_year"]
        ),
        demand=demand_file,
        costs=technology_data("../technology-data/outputs/costs_{year}.csv"),
        wacc="data/wacc.csv",
        network="resources/networks_ip/{scenario}/{year}/{esc}/{from}-{to}/network.nc",
        additional_components="resources/additional_components.pkl",
    output:
        network=(
            "resources/networks_ip_as/{scenario}/{year}/{esc}/{from}-{to}/network.nc"
        ),
        lcoes="resources/networks_ip_as/{scenario}/{year}/{esc}/{from}-{to}/lcoes.csv",
    threads: 1
    params:
        scenario=lambda w: get_scenario(w["scenario"]),
    log:
        python="logs/{scenario}/{year}/{esc}/{from}-{to}/attach_supply.log",
        notebook="logs/{scenario}/{year}/{esc}/{from}-{to}/attach_supply.ipynb",
    notebook:
        "../actions/attach_supply.py.ipynb"

# SPDX-FileCopyrightText: 2020-2022 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later


rule download_technology_data:
    input:
        HTTP.remote(
            expand(
                "raw.githubusercontent.com/PyPSA/technology-data/{version}/outputs/costs_{year}.csv",
                version=config["technology_data"],
                allow_missing=True,
            ),
            keep_local=True,
        ),
    output:
        "data/technology-data/outputs/costs_{year}.csv",
    run:
        move(input[0], output[0])


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
        costs="data/technology-data/outputs/costs_{year}.csv",
        wacc="data/wacc.csv",
        distances="data/distances.csv",
        shipping_properties="data/shipping.csv",
        network="escs/{esc}",
        additional_components="resources/additional_components.pkl",
    output:
        network="resources/networks/{scenario}/{year}/{esc}/{from}-{to}/network.nc",
    group: "esc"
    params:
        scenario=lambda w: get_scenario(w["scenario"]),
        era_year=config["renewables"]["era5_year"],
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
    group: "esc"
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


rule attach_supply:
    message:
        "Attaching RES supply ('as') to network."
    input:
        supply="resources/supply_{from}.nc",
        demand="data/overwrite/demand.csv",
        costs="data/technology-data/outputs/costs_{year}.csv",
        wacc="data/wacc.csv",
        efficiencies="data/efficiencies.csv",
        network="resources/networks_ip/{scenario}/{year}/{esc}/{from}-{to}/network.nc",
        additional_components="resources/additional_components.pkl",
    output:
        network=(
            "resources/networks_ip_as/{scenario}/{year}/{esc}/{from}-{to}/network.nc"
        ),
        lcoes="results/{scenario}/{year}/{esc}/{from}-{to}/lcoes.csv",
    group: "esc"
    threads: 1
    params:
        scenario=lambda w: get_scenario(w["scenario"]),
        renewable_details=config["renewables"],
    log:
        python="logs/{scenario}/{year}/{esc}/{from}-{to}/attach_supply.log",
        notebook="logs/{scenario}/{year}/{esc}/{from}-{to}/attach_supply.ipynb",
    notebook:
        "../actions/attach_supply.py.ipynb"

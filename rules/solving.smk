# SPDX-FileCopyrightText: 2020-2021 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later

rule solve_scenario:
    input:
        expand("results/{scenario}/{year}/{esc}/{exporter}-{importer}/network.nc",
                scenario=SCENARIO,
                esc=ESCS,
                exporter=EXPORTERS,
                importer=IMPORTERS,
                year=YEARS
                ),
        inputs=f"results/{SCENARIO}/inputs.tar"
        
rule solve_network:
    input:
        network="resources/networks_supplied/{scenario}/{year}/{esc}/{from}-{to}/network.nc",
        additional_components="resources/additional_components.pkl"
    output:
        network="results/{scenario}/{year}/{esc}/{from}-{to}/network.nc"
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 5000
    log:
        python="logs/{scenario}/{year}/{esc}/{from}-{to}/solve_network.log",
        notebook="logs/{scenario}/{year}/{esc}/{from}-{to}/solve_network.ipynb"
    notebook:
        "../actions/solve_network.py.ipynb"

rule backup_scenario:
    input:
        config="./config.yaml",
        data="./data/",
        costs="../technology-data/outputs/"
    output:
        tarchive=f"results/{SCENARIO}/inputs.tar",
    threads: 1
    script:
        "../actions/backup_run.py"

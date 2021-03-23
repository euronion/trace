# SPDX-FileCopyrightText: 2020-2021 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later

rule solve_all:
    input:
        expand("results/"+SCENARIO_FOLDER+"/{esc}/{exporter}-{importer}/network.nc", esc=ESCS, exporter=EXPORTERS, importer=IMPORTERS),
        config="results/"+SCENARIO_FOLDER+"/config.yaml"
        
rule solve_network:
    input:
        network="resources/networks_supplied/"+SCENARIO_FOLDER+"/{esc}/{from}-{to}/network.nc",
        additional_components="resources/additional_components.pkl"
    output:
        network="results/"+SCENARIO_FOLDER+"/{esc}/{from}-{to}/network.nc"
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 5000
    log:
        python="logs/"+SCENARIO_FOLDER+"/solve_network/{esc}/{from}-{to}.log",
        notebook="logs/"+SCENARIO_FOLDER+"/solve_network/{esc}/{from}-{to}.ipynb"
    notebook:
        "../actions/solve_network.py.ipynb"

rule backup_run:
    input:
        config="config.yaml",
        data="data/",
        costs=f"../technology-data/outputs/costs_{config['scenario']['year']}.csv"
    output:
        config="results/"+SCENARIO_FOLDER+"/config.yaml",
        data="results/"+SCENARIO_FOLDER+"/data.tar",
        costs=f"results/{SCENARIO_FOLDER}/costs_{config['scenario']['year']}.csv"
    threads: 1
    script:
        "../actions/backup_run.py"

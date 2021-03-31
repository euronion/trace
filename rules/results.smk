# SPDX-FileCopyrightText: 2020-2021 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later

def all_result_files(wildcards):

    return expand("results/{scenario}/{year}/{esc}/{exporter}-{importer}/results.csv",
                scenario=wildcards.scenario,
                year=config['scenarios'][wildcards.scenario]['YEARS'],
                esc=config['scenarios'][wildcards.scenario]['ESCS'],
                exporter=config['scenarios'][wildcards.scenario]['EXPORTERS'],
                importer=config['scenarios'][wildcards.scenario]['IMPORTERS'],
                )

rule combine_scenario_results:
    input:
        results=all_result_files,
        inputs="results/{scenario}/inputs.tar"
    output:
        results="results/{scenario}/results.csv"
    threads: 1
    log:
        python="logs/{scenario}/combine_results.log",
        notebook="logs/{scenario}/combine_results.ipynb"
    notebook:
        "../actions/combine_results.py.ipynb"
        
rule extract_result:
    input:
        network="results/{scenario}/{year}/{esc}/{exporter}-{importer}/network.nc"
    output:
        results="results/{scenario}/{year}/{esc}/{exporter}-{importer}/results.csv"
    threads: 1
    log:
        python="logs/{scenario}/{year}/{esc}/{exporter}-{importer}/extract_result.log",
        notebook="logs/{scenario}/{year}/{esc}/{exporter}-{importer}/extract_result.ipynb",
    notebook:
        "../actions/extract_result.py.ipynb"

# SPDX-FileCopyrightText: 2020-2021 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later

rule combine_scenario_results:
    input:
        expand("results/{scenario}/{year}/{esc}/{exporter}-{importer}/results.csv",
                esc=ESCS,
                exporter=EXPORTERS,
                importer=IMPORTERS,
                year=YEARS,
                allow_missing=True
                ),
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

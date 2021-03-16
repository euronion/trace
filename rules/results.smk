# SPDX-FileCopyrightText: 2020-2021 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later

rule combine_results:
    input:
        expand("results/{year}_{wacc}/{esc}/{exporter}-{importer}/results.csv", year=YEARS, wacc=WACCS, esc=ESCS, exporter=EXPORTERS, importer=IMPORTERS)
    output:
        results="results/results.csv"
    threads: 1
    log:
        python="logs/combine_results.log",
        notebook="logs/combine_results.ipynb"
    notebook:
        "actions/combine_results.py.ipynb"
        
rule extract_result:
    input:
        network="results/{year}_{wacc}/{esc}/{from}-{to}/network.nc"
    output:
        results="results/{year}_{wacc}/{esc}/{from}-{to}/results.csv"
    threads: 1
    log:
        python="logs/{year}_{wacc}/extract_result/{esc}/{from}-{to}.log",
        notebook="logs/{year}_{wacc}/extract_result/{esc}/{from}-{to}.ipynb"
    notebook:
        "actions/extract_result.py.ipynb"

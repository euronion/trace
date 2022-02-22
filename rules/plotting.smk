# SPDX-FileCopyrightText: 2020-2022 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later

rule plot_supply_curves:
    input:
        lcoes="results/{scenario}/{year}/{esc}/{exporter}-{importer}/lcoes.csv",
    output:
        fig=multiext("results/{scenario}/{year}/{esc}/{exporter}-{importer}/supply-curves", ".png", ".html")
    threads: 1
    params:
        scenario=lambda w: get_scenario(w["scenario"]),
    log:
        python="logs/{scenario}/{year}/{esc}/{exporter}-{importer}/plot_supply_curves.log",
        notebook=(
            "logs/{scenario}/{year}/{esc}/{exporter}-{importer}/plot_supply_curves.ipynb"
        ),
    notebook:
        "../actions/plot_supply_curves.py.ipynb"

rule plot_generation:
    message: "Plotting generation PoV: RES CFs, capacities, cost and electricity generation."
    input:
        network="results/{scenario}/{year}/{esc}/{exporter}-{importer}/network.nc",
    output:
        fig=multiext("results/{scenario}/{year}/{esc}/{exporter}-{importer}/generation", ".png", ".html")
    threads: 1
    params:
        scenario=lambda w: get_scenario(w["scenario"]),
    log:
        python="logs/{scenario}/{year}/{esc}/{exporter}-{importer}/plot_generation.log",
        notebook=(
            "logs/{scenario}/{year}/{esc}/{exporter}-{importer}/plot_generation.ipynb"
        ),
    notebook:
        "../actions/plot_generation.py.ipynb"

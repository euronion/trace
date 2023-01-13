# SPDX-FileCopyrightText: 2023 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later


rule plot_LCoEs_by_flexibility:
    input:
        results="results/results.csv",
    output:
        pdf="figures/LCoE_by_flexibility/{esc}_{exporter}.pdf",
        png="figures/LCoE_by_flexibility/{esc}_{exporter}.png",
    log:
        "logs/plot_LCoEs_by_flexibility/{esc}_{exporter}.log",
    notebook:
        "../actions/paper_plots/plot_LCoEs_by_flexibility.py.ipynb"


rule plot_LCoEs_per_exporter_by_flexibility:
    input:
        results="results/results.csv",
    output:
        pdf="figures/LCoEs_per_exporter_by_flexibility/{exporter}.pdf",
        png="figures/LCoEs_per_exporter_by_flexibility/{exporter}.png",
    log:
        "logs/plot_LCoEs_per_exporter_by_flexibility/{exporter}.log",
    notebook:
        "../actions/paper_plots/plot_LCoEs_per_exporter_by_flexibility.py.ipynb"


rule plot_storage_link_utilisation:
    input:
        results="results/results.csv",
    output:
        pdf="figures/storage_link_utilisation/{esc}_{exporter}.pdf",
        png="figures/storage_link_utilisation/{esc}_{exporter}.png",
    log:
        "logs/plot_storage_link_utilisation/{esc}_{exporter}.log",
    notebook:
        "../actions/paper_plots/plot_storage_link_utilisation.py.ipynb"


rule plot_RES_shares:
    input:
        results="results/results.csv",
    output:
        pdf="figures/RES_shares/{esc}_{exporter}.pdf",
        png="figures/RES_shares/{esc}_{exporter}.png",
    log:
        python="logs/plot_RES_shares/{esc}_{exporter}.log",
    notebook:
        "../actions/paper_plots/plot_RES_shares.py.ipynb"


rule plot_cost_compositions:
    input:
        results="results/results.csv",
    output:
        pdf="figures/cost_composition/{esc}_{exporter}.pdf",
        png="figures/cost_composition/{esc}_{exporter}.png",
    log:
        python="logs/cost_composition/{esc}_{exporter}.log",
    notebook:
        "../actions/paper_plots/plot_cost_compositions.py.ipynb"


rule plot_import_route:
    input:
        gebco="resources/gebco/GEBCO_2021.nc",
    output:
        pdf="figures/import_route.pdf",
        png="figures/import_route.png",
    log:
        "logs/plot_import_route.log",
    notebook:
        "../actions/paper_plots/plot_import_route.py.ipynb"


rule plot_all_paper_figures:
    default_target: True
    input:
        expand(
            rules.plot_LCoEs_by_flexibility.output.pdf,
            esc=["hvdc-to-elec", "pipeline-h2-to-elec"],
            exporter=["MA", "TN"],
        ),
        expand(
            rules.plot_LCoEs_per_exporter_by_flexibility.output.pdf,
            exporter=["MA", "TN"],
        ),
        expand(
            rules.plot_storage_link_utilisation.output.pdf,
            esc=["hvdc-to-elec", "pipeline-h2-to-elec"],
            exporter=["MA", "TN"],
        ),
        expand(
            rules.plot_RES_shares.output.pdf,
            esc=["hvdc-to-elec", "pipeline-h2-to-elec"],
            exporter=["MA", "TN"],
        ),
        rules.plot_import_route.output.pdf,

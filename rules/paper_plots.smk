# SPDX-FileCopyrightText: 2023 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later


rule plot_figA:
    input:
        results="results/results.csv",
    output:
        pdf="figures/figA/figA.pdf",
        png="figures/figA/figA.png",
    log:
        "logs/figA/figA.log",
    notebook:
        "../actions/paper_plots/plot_figA.py.ipynb"


rule plot_figB:
    input:
        results="results/results.csv",
    output:
        pdf="figures/figB/{exporter}_{csp}.pdf",
        png="figures/figB/{exporter}_{csp}.png",
    log:
        "logs/figB/{exporter}_{csp}.log",
    wildcard_constraints:
        csp="with_csp|without_csp",
    notebook:
        "../actions/paper_plots/plot_figB.py.ipynb"


rule plot_figC:
    input:
        results="results/results.csv",
    output:
        pdf="figures/figC/{exporter}_{csp}.pdf",
        png="figures/figC/{exporter}_{csp}.png",
    log:
        "logs/figC/{exporter}_{csp}.log",
    wildcard_constraints:
        csp="with_csp|without_csp",
    notebook:
        "../actions/paper_plots/plot_figC.py.ipynb"


rule plot_figD:
    input:
        results="results/results.csv",
    output:
        pdf="figures/figD/{esc}_{exporter}.pdf",
        png="figures/figD/{esc}_{exporter}.png",
    log:
        "logs/figD/{esc}_{exporter}.log",
    notebook:
        "../actions/paper_plots/plot_figD.py.ipynb"


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
        pdf="figures/storage_link_utilisation/{exporter}_{csp}.pdf",
        png="figures/storage_link_utilisation/{exporter}_{csp}.png",
    log:
        "logs/plot_storage_link_utilisation/{exporter}_{csp}.log",
    wildcard_constraints:
        csp="with_csp|without_csp",
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
            exporter=["MA", "TN"],
            csp=["with_csp", "without_csp"],
        ),
        expand(
            rules.plot_RES_shares.output.pdf,
            esc=["hvdc-to-elec", "pipeline-h2-to-elec"],
            exporter=["MA", "TN"],
        ),
        rules.plot_figA.output.pdf,
        expand(
            rules.plot_figB.output.pdf,
            exporter=["MA", "TN"],
            csp=["with_csp", "without_csp"],
        ),
        expand(
            rules.plot_figC.output.pdf,
            exporter=["MA", "TN"],
            csp=["with_csp", "without_csp"],
        ),
        expand(
            rules.plot_figD.output.pdf,
            esc=["hvdc-to-elec", "pipeline-h2-to-elec"],
            exporter=["MA", "TN"],
        ),
        rules.plot_import_route.output.pdf,

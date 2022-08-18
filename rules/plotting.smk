# SPDX-FileCopyrightText: 2020-2022 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later


rule plot_LCoEs_single_year:
    input:
        results=rules.all_scenario_results.output,
    output:
        figures=multiext(
            "figures/paper-01/LCoE_all-ESC-EXP_{year}_{scenario}", ".pdf", ".png"
        ),
    notebook:
        "../actions/plotting/LCoE_all-ESC-EXP.py.ipynb"


rule plot_2030_to_2050_LCoEs:
    input:
        results=rules.all_scenario_results.output,
    output:
        figures=multiext(
            "figures/paper-01/LCoE_all-ESC-EXP_2030-2050_{scenario}", ".pdf", ".png"
        ),
    notebook:
        "../actions/plotting/LCoE_all-ESC-EXP_2030-2050.py.ipynb"


rule plot_2030_to_2050_LCoHs:
    input:
        results=rules.all_scenario_results.output,
    output:
        figures=multiext(
            "figures/paper-01/LCoH_all-ESC-EXP_2030-2050_{scenario}", ".pdf", ".png"
        ),
    notebook:
        "../actions/plotting/LCoH_all-ESC-EXP_2030-2050.py.ipynb"


rule plot_selected_ESCs_cost_compositions:
    input:
        results=rules.all_scenario_results.output,
    output:
        figures=multiext(
            "figures/paper-01/cost-compositions_selected_ESCs_{scenario}_{year}",
            ".pdf",
            ".png",
        ),
    notebook:
        "../actions/plotting/cost-compositions_selected_ESCs.py.ipynb"


rule plot_sensitivities:
    input:
        results=rules.all_scenario_results.output,
    output:
        figures=multiext(
            "figures/paper-01/sensitivities_{year}_{esc}_{exporter}-{importer}",
            ".pdf",
            ".png",
        ),
    notebook:
        "../actions/plotting/sensitivities_selected_esc.py.ipynb"


rule plot_electricity_supply_curves:
    input:
        supply_curves=expand(
            "resources/networks_ip_as/{scenario}/{year}/hvdc-to-h2/{exp}-DE/lcoes.csv",
            exp=["AR", "AU", "DE", "DK", "EG", "ES", "MA", "SA"],
            allow_missing=True,
        ),
        domestic_demand="resources/demand_annual_TRACES_2013.csv",
    output:
        figures=multiext(
            "figures/paper-01/electricity_supply-curves_{year}_{scenario}",
            ".pdf",
            ".png",
        ),
    notebook:
        "../actions/plotting/electricity_supply-curves.py.ipynb"


rule plot_electricity_generation_shares:
    input:
        results=rules.all_scenario_results.output,
    output:
        figures=multiext(
            "figures/paper-01/res_generation_shares_all-ESC-EXP_{year}_{scenario}",
            ".pdf",
            ".png",
        ),
    notebook:
        "../actions/plotting/res_generation_shares_all-ESC-EXP.py.ipynb"


rule plot_ESF_vs_LCoEs_vs_curtailment:
    input:
        results=rules.all_scenario_results.output,
    output:
        figures=multiext(
            "figures/paper-01/ESF_vs_LCoEs_vs_curtailment_{year}_{scenario}",
            ".pdf",
            ".png",
        ),
    notebook:
        "../actions/plotting/ESF_vs_LCoEs_vs_curtailment.py.ipynb"


rule plot_all_paper_figures:
    input:
        expand(
            "figures/paper-01/LCoE_all-ESC-EXP_{year}_{scenario}.pdf",
            year=[2030, 2040, 2050],
            scenario=["default", "lowhomogeneous"],
        ),
        expand(
            "figures/paper-01/LCoE_all-ESC-EXP_2030-2050_{scenario}.pdf",
            scenario=["default", "lowhomogeneous"],
        ),
        expand(
            "figures/paper-01/LCoH_all-ESC-EXP_2030-2050_{scenario}.pdf",
            scenario=["default", "lowhomogeneous"],
        ),
        # Sensivity cases as included in scenarios.csv
        expand(
            "figures/paper-01/sensitivities_2030_{esc}_ES-DE.pdf",
            esc=["pipeline-h2", "shipping-meoh"],
        ),
        # Cost compositions
        expand(
            "figures/paper-01/cost-compositions_selected_ESCs_{scenario}_{year}.pdf",
            scenario=["default", "lowhomogeneous"],
            year=[2030, 2040, 2050],
        ),
        expand(
            "figures/paper-01/electricity_supply-curves_{year}_{scenario}.pdf",
            year=[2030, 2040, 2050],
            scenario=["default", "lowhomogeneous"],
        ),
        expand(
            "figures/paper-01/res_generation_shares_all-ESC-EXP_{year}_{scenario}.pdf",
            year=[2030, 2040, 2050],
            scenario=["default", "lowhomogeneous"],
        ),
        expand(
            "figures/paper-01/ESF_vs_LCoEs_vs_curtailment_{year}_{scenario}.pdf",
            year=[2030, 2040, 2050],
            scenario=["default", "lowhomogeneous"],
        ),

# SPDX-FileCopyrightText: 2020-2022 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later


rule plot_2030_LCoEs:
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
            "figures/paper-01/cost-compositions_selected_ESCs_{scenario}",
            ".pdf",
            ".png",
        ),
    notebook:
        "../actions/plotting/cost-compositions_selected_ESCs.py.ipynb"  #TODO


rule plot_sensitivities:
    input:
        results=rules.all_scenario_results.output,
    output:
        figures=multiext(
            "figures/paper-01/sensitivities_{year}_{esc}_{exporter}-{importer}", ".pdf", ".png"
        ),
    notebook:
        "../actions/plotting/sensitivities_selected_esc.py.ipynb"


rule plot_all_paper_figures:
    input:
        rules.plot_2030_LCoEs.output,
        rules.plot_2030_to_2050_LCoEs.output,
        rules.plot_2030_to_2050_LCoHs.output,
        rules.plot_selected_ESCs_cost_compositions.output,
        # Sensivity cases as included in scenarios.csv
        "figures/paper-01/sensitivities_2030_pipeline-h2_ES-DE.pdf",
        "figures/paper-01/sensitivities_2030_shipping-meoh_ES-DE.pdf",

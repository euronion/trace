# SPDX-FileCopyrightText: 2020-2022 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later

# Use paramspace to evaluate which scenarios to run
scenarios = Paramspace(pd.read_csv("scenarios/default.csv", dtype=str))


# Custom pattern for formatting Paramspace, as snakemake
# does currently not allow for patterns without the wildcard_name included
# see: https://stackoverflow.com/questions/71293563/custom-patterns-with-snakemakes-paramspace/71296522#71296522
def custom_instance_pattern(ps):
    pattern = "{scenario}/{year}/{esc}/{exporter}-{importer}"
    instance_patterns = [
        pattern.format(**dict(i for i in row.items())) for _, row in ps.iterrows()
    ]
    return instance_patterns


rule all_scenario_results:
    input:
        results=expand(
            "results/{instances}/results.csv",
            instances=custom_instance_pattern(scenarios),
        ),
    output:
        results="results/results.csv",
    log:
        python="logs/combine_results.log",
        notebook="logs/combine_results.ipynb",
    threads: 1
    notebook:
        "../actions/combine_results.py.ipynb"


rule extract_result:
    input:
        network="results/{scenario}/{year}/{esc}/{exporter}-{importer}/network.nc",
    output:
        results="results/{scenario}/{year}/{esc}/{exporter}-{importer}/results.csv",
    log:
        python="logs/{scenario}/{year}/{esc}/{exporter}-{importer}/extract_result.log",
        notebook=(
            "logs/{scenario}/{year}/{esc}/{exporter}-{importer}/extract_result.ipynb"
        ),
    threads: 1
    params:
        scenario=lambda w: get_scenario(w["scenario"]),
    notebook:
        "../actions/extract_result.py.ipynb"

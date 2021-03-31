# SPDX-FileCopyrightText: 2020-2021 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later

configfile: "config.yaml"

wildcard_constraints:
    year="\d+",
    scenario="[-\w]+"
    
#if 'scenario' not in globals():
#    raise KeyError("Need to define 'scenario' as parameter "
#                   "using the --config scenario=<...> flag for snakemake.")
                       
#SCENARIO=scenario
#ESCS=config["scenarios"][SCENARIO]["ESCS"]
#EXPORTERS=config["scenarios"][SCENARIO]["EXPORTERS"]
#IMPORTERS=config["scenarios"][SCENARIO]["IMPORTERS"]
#YEARS=config["scenarios"][SCENARIO]["YEARS"]

subworkflow technology_data:
    workdir:
        "../technology-data"
    snakefile:
        "../technology-data/Snakefile"
    configfile:
        "../technology-data/config.yaml"

include: "rules/gegis.smk"
include: "rules/esc_construction.smk"
include: "rules/solving.smk"
include: "rules/results.smk"

# SPDX-FileCopyrightText: 2020-2021 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later


configfile: "config.yaml"


wildcard_constraints:
    year="\d+",
    scenario="[-\w]+",


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

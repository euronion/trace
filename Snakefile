# SPDX-FileCopyrightText: 2020-2021 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later

configfile: "config.yaml"

ESCS=["hvdc","pipeline-h2","pipeline-ch4","shipping-lh2","shipping-lch4","shipping-meoh","shipping-lnh3","shipping-lohc"]
EXPORTERS=["AU","AR","ES","EG","MA","SA","DK","DE"]
IMPORTERS=["DE"]
YEARS=[2030,2040,2050]
WACCS=["homogeneous", "lowhomogeneous"]

SCENARIO_FOLDER = f"{config['scenario']['year']}_{config['scenario']['wacc']}"

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


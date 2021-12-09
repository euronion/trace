# SPDX-FileCopyrightText: 2020-2021 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later

from pathlib import Path
from snakemake.utils import update_config
from snakemake.io import load_configfile


# Specify config file
configfile: "config/config.cern_link.yaml"


# Default configs - do not change
default_configfile = "config/config.default.yaml"

# Load default config and overwrite with specific config
specific_config = config.copy()
config = load_configfile(Path(default_configfile))
update_config(config, specific_config)


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


def get_scenario(scenario_name):
    s = config["scenarios"]["default"].copy()
    update_config(s, config["scenarios"][scenario_name])
    return s


include: "rules/gegis.smk"
include: "rules/esc_construction.smk"
include: "rules/solving.smk"
include: "rules/results.smk"

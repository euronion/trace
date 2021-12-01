# SPDX-FileCopyrightText: 2020-2021 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later

from pathlib import Path
from snakemake.utils import update_config
from snakemake.io import load_configfile

# Specify config file
configfile: "config/config.initial_paper.yaml"

# Default configs - do not change
default_configfile: "config/config.default.yaml"

# Load default config and overwrite with specific config
config = load_configfile(default_configfile)
update_config(config, load_configfile(configfile))

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

# SPDX-FileCopyrightText: 2020-2022 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later

from pathlib import Path
from copy import deepcopy
import pandas as pd
from snakemake.utils import Paramspace
from snakemake.utils import update_config
from shutil import move
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from snakemake.io import load_configfile

HTTP = HTTPRemoteProvider()


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


def get_scenario(scenario_name):
    s = deepcopy(config["scenarios"]["default"])
    update_config(s, config["scenarios"][scenario_name])
    return s


include: "rules/renewables.smk"
include: "rules/esc_construction.smk"
include: "rules/solving.smk"
include: "rules/results.smk"
include: "rules/plotting.smk"

#!/usr/bin/env python

# SPDX-FileCopyrightText: 2020-2021 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later

# coding: utf-8

from pathlib import Path
import shutil
import tarfile
import tempfile
import yaml
import json

if __name__ == "__main__":

    # Write the complete config into a temporary config file
    temp_dir = tempfile.TemporaryDirectory()
    config_path = Path(temp_dir.name) / "config.yaml"

    # Recursively convert OrderedDicts to Dicts (OrderedDicts are not supported by yaml.safe_dump)
    config = json.loads(json.dumps(snakemake.config))
    with open(config_path, "w") as config_file:
        yaml.safe_dump(config, config_file, default_flow_style=False, sort_keys=False)

    with tarfile.open(snakemake.output["tarchive"], "w") as tararchive:
        for k, src_path in snakemake.input.items():
            tararchive.add(src_path, recursive=True)

        # Add temp config file
        tararchive.add(config_path, arcname="config.complete.yaml")

    # Remove temporary directory used for config file dumping
    temp_dir.cleanup()

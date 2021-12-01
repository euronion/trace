#!/usr/bin/env python

# SPDX-FileCopyrightText: 2020-2021 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later

# coding: utf-8

from pathlib import Path
import shutil
import tarfile

if __name__ == "__main__":

    with tarfile.open(snakemake.output["tarchive"], "w") as tararchive:
        for k, src_path in snakemake.input.items():
            tararchive.add(src_path, recursive=True)

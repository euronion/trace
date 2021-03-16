#!/usr/bin/env python

# SPDX-FileCopyrightText: 2020-2021 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later

# coding: utf-8

from pathlib import Path
import shutil
import tarfile

if __name__ == '__main__':
    
    for k, src_path in snakemake.input.items():
        src_path = Path(src_path)
        tar_path = Path(snakemake.output[k])
        
        if src_path.is_dir():
            with tarfile.open(tar_path, 'w') as archive:
                archive.add(src_path, recursive=True)
        else:
            shutil.copy(src_path, snakemake.output[k])
# SPDX-FileCopyrightText: 2022 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later

# - Rules for determining renewables potentials and time-series - #

from pathlib import Path

# For downloading files using Snakemake
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

# rule build_cutout:


# Downloading GADM database for country/region shapes which are used to define
# export region extents.
# Website: https://gadm.org/
rule download_gadm:
    input:
        HTTP.remote(
            "biogeo.ucdavis.edu/data/gadm3.6/gadm36_gpkg.zip",
            static=True,
        ),
    output:
        "resources/GADM/gadm36.gpkg",
    run:
        output_folder = Path(output[0]).parent
        shell("unzip {input} -d {output_folder}")


# rule download_landcover:
# rule download_protected_areas:
# rule download_population_map:

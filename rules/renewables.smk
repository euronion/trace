# SPDX-FileCopyrightText: 2022 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later

# - Rules for determining renewables potentials and time-series - #

from pathlib import Path

# For downloading files using Snakemake
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

rule download_eez:
    message:
        """Trying to download EEZ files. If this fail, manually download from https://www.marineregions.org/download_file.php?name=World_EEZ_v11_20191118_gpkg.zip and extract them into the 'resources/' folder."""
    output:
        zip="resources/World_EEZ_v11_20191118_gpkg.zip",
        gpkg="resources/World_EEZ_v11_20191118_gpkg/eez_v11.gpkg",
    run:
        shell("curl  -X POST --data 'name=Name&organisation=Organisation&email=e.mail%40inter.net&country=Germany&user_category=academia&purpose_category=Research&agree=1' 'https://www.marineregions.org/download_file.php?name=World_EEZ_v11_20191118_gpkg.zip' --output '{output.zip}'")
        output_folder = Path(output["zip"]).parent
        shell("unzip {output.zip} -d {output_folder}")

rule build_region_shape:
    message:
        "Creating region definition for: {wildcards.region}."
    input:
        gadm="resources/gadm/gadm36_levels.gpkg",
        eez="resources/World_EEZ_v11_20191118_gpkg/eez_v11.gpkg",
    output:
        gpkg="resources/regions/{region}.gpkg"
    params:
        region_members=lambda w: config["regions"][w["region"]]
    notebook:
        "../actions/build_region_shape.py.ipynb"

rule build_cutout:
    message:
        "Creating cutout: {output}."
    input:
        gadm="resources/gadm/gadm36.gpkg"
    output:
        cutout="resources/cutouts/{region}.nc"
    notebook:
        "../actions/build_cutout.py.ipynb"


# Downloading GADM database for country/region shapes which are used to define
# export region extents.
# Website: https://gadm.org/
rule download_gadm:
    input:
        HTTP.remote(
            "biogeo.ucdavis.edu/data/gadm3.6/gadm36_levels_gpkg.zip",
            static=True,
        ),
    output:
        "resources/gadm/gadm36_levels.gpkg",
    run:
        output_folder = Path(output[0]).parent
        shell("unzip {input} -d {output_folder}")

# rule download_landcover:
# rule download_protected_areas:
# rule download_population_map:

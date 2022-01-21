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

# Downloading Copernicus Global Land Cover for land cover and land use:
# Website: https://land.copernicus.eu/global/products/lc
rule download_land_cover:
    input:
        HTTP.remote(
            "zenodo.org/record/3939050/files/PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif?download=1",
            static=True,
        ),
    output:
        "resources/Copernicus_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif",
    shell:
        "mv {input} {output}"

# Downloading bathymetry data (GEBCO)
# Website: https://www.gebco.net/data_and_products/gridded_bathymetry_data/#global
rule download_gebco:
    input:
        HTTP.remote(
            "www.bodc.ac.uk/data/open_download/gebco/gebco_2021_tid/zip",
            additional_request_string="/",
            static=True,
        ),
    output:
        zip="resources/gebco_2021_tid.zip",
        gebco="resources/gebco/GEBCO_2021_TID.nc",
    run:
        shell("mv {input} {output.zip}")
        output_folder = Path(output["gebco"]).parent
        shell("unzip {output.zip} -d {output_folder}")

# Downloading Marine protected area database from WDPA
# Website: https://www.protectedplanet.net/en/thematic-areas/marine-protected-areas
rule download_wdpa_marine:
    input:
        HTTP.remote(
            "d1gam3xoknrgr2.cloudfront.net/current/WDPA_WDOECM_Jan2022_Public_marine_shp.zip",
            static=True,
        ),
    output:
        zip="resources/WDPA_WDOECM_Jan2022_marine_shp.zip",
        folder=directory("resources/wdpa_marine/"),
    run:
        shell("mv {input} {output.zip}")
        shell("unzip {output.zip} -d {output.folder}")

# Downloading Marine protected area database from WDPA
# Website: https://www.protectedplanet.net/en/thematic-areas/wdpa
rule download_wdpa:
    input:
        HTTP.remote(
            "d1gam3xoknrgr2.cloudfront.net/current/WDPA_Jan2022_Public_shp.zip",
            static=True,
        ),
    output:
        zip="resources/WDPA_Jan2022_shp.zip",
        folder=directory("resources/wdpa/"),
    run:
        shell("mv {input} {output.zip}")
        shell("unzip {output.zip} -d {output.folder}")

rule build_region_shape:
    message:
        "Creating region definition (off and onshore) for: {wildcards.region}."
    input:
        gadm="resources/gadm/gadm36_levels.gpkg",
        eez="resources/World_EEZ_v11_20191118_gpkg/eez_v11.gpkg",
    output:
        gpkg="resources/regions/{region}.gpkg",
    params:
        region_members=lambda w: config["regions"][w["region"]],
        offshore_proximity=lambda w: config["regions"]["offshore_proximity"],
    notebook:
        "../actions/build_region_shape.py.ipynb"

rule build_region_cutout:
    message:
        "Downloading cutout for: {wildcards.region}."
    input:
        gpkg="resources/regions/{region}.gpkg",
    output:
        cutout="resources/cutouts/{region}.nc"
    params:
        era5_year=config["renewables"]["era5_year"]
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

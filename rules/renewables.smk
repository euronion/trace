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
        shell(
            "curl  -X POST --data 'name=Name&organisation=Organisation&email=e.mail%40inter.net&country=Germany&user_category=academia&purpose_category=Research&agree=1' 'https://www.marineregions.org/download_file.php?name=World_EEZ_v11_20191118_gpkg.zip' --output '{output.zip}'"
        )
        output_folder = Path(output["zip"]).parent
        shell("unzip {output.zip} -d {output_folder}")


# Downloading Copernicus Global Land Cover for land cover and land use:
# Website: https://land.copernicus.eu/global/products/lc
rule download_land_cover:
    input:
        HTTP.remote(
            "zenodo.org/record/3939050/files/PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif",
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
            "www.bodc.ac.uk/data/open_download/gebco/gebco_2021/zip",
            additional_request_string="/",
            static=True,
        ),
    output:
        zip="resources/gebco_2021.zip",
        gebco="resources/gebco/GEBCO_2021.nc",
    run:
        shell("mv {input} {output.zip}")
        output_folder = Path(output["gebco"]).parent
        shell("unzip {output.zip} -d {output_folder}")


# Transform the GEBCO dataset into a second dataset
# with georeferenced slope values (in percent) for all locations
# Projection of the output dataset is Global Mollweide (ESRI:54009)
rule calcualte_gebco_slope:
    input:
        gebco="resources/gebco/GEBCO_2021.nc",
    output:
        mollweide=temp("resources/gebco/GEBCO_2021_mollweide.nc"),
        slope_inflated=temp("resources/gebco/GEBCO_2021_mollweide_inflated.nc"),
        slope="resources/gebco/GEBCO_2021_slope.nc",
    run:
        shell(
            "gdalwarp -multi -wo NUM_THREADS={threads} -of netCDF -co FORMAT=NC4 -s_srs 'EPSG:4326' -t_srs 'ESRI:54009' {input.gebco} {output.mollweide}"
        )
        shell(
            "gdaldem slope -p -of netCDF -co FORMAT=NC4 {output.mollweide} {output.slope_inflated}"
        )
        shell(
            "gdal_translate -of netCDF -co FORMAT=NC4 -co COMPRESS=DEFLATE -co ZLEVEL=1 {output.slope_inflated} {output.slope}"
        )


# Downloading Marine protected area database from WDPA
# extract the main zip and then merge the contained 3 zipped shapefiles
# Website: https://www.protectedplanet.net/en/thematic-areas/marine-protected-areas
rule download_wdpa_marine:
    input:
        HTTP.remote(
            "d1gam3xoknrgr2.cloudfront.net/current/WDPA_WDOECM_Feb2022_Public_marine_shp.zip",
            static=True,
        ),
    output:
        zip="resources/WDPA_WDOECM_Feb2022_marine.zip",
        folder=directory("resources/WDPA_WDOECM_Feb2022_marine"),
        gpkg="resources/WDPA_WDOECM_Feb2022_marine.gpkg",
    run:
        shell("mv {input} {output.zip}")
        shell("unzip {output.zip} -d {output.folder}")
        for i in range(3):
            # vsizip is special driver for directly working with zipped shapefiles in ogr2ogr
            layer_path = (
                f"/vsizip/{output.folder}/WDPA_WDOECM_Feb2022_Public_marine_shp_{i}.zip"
            )
            print(f"Adding layer {i+1} of 3 to combined output file.")
            shell("ogr2ogr -f gpkg -update -append {output.gpkg} {layer_path}")


# Downloading protected area database from WDPA
# extract the main zip and then merge the contained 3 zipped shapefiles
# Website: https://www.protectedplanet.net/en/thematic-areas/wdpa
rule download_wdpa:
    input:
        HTTP.remote(
            "d1gam3xoknrgr2.cloudfront.net/current/WDPA_Feb2022_Public_shp.zip",
            static=True,
        ),
    output:
        zip="resources/WDPA_Feb2022_shp.zip",
        folder=directory("resources/WDPA_Feb2022"),
        gpkg="resources/WDPA_Feb2022.gpkg",
    run:
        shell("mv {input} {output.zip}")
        shell("unzip {output.zip} -d {output.folder}")
        for i in range(3):
            # vsizip is special driver for directly working with zipped shapefiles in ogr2ogr
            layer_path = f"/vsizip/{output.folder}/WDPA_Feb2022_Public_shp_{i}.zip"
            print(f"Adding layer {i+1} of 3 to combined output file.")
            shell("ogr2ogr -f gpkg -update -append {output.gpkg} {layer_path}")


# Downloading global shipping traffic density map
# with ship movement between 2015-2020
# Specific dataset: "Global Ship Density"
# Website: https://datacatalog.worldbank.org/search/dataset/0037580/Global-Shipping-Traffic-Density
rule download_shipping_density:
    output:
        zip="resources/shipdensity_global.zip",
        shipdensity="resources/shipdensity/shipdensity_global.tif",
    run:
        shell(
            "curl 'https://datacatalogapi.worldbank.org/ddhxext/ResourceDownload?resource_unique_id=DR0045406&version_id=2022-06-28T08:26:08.9439370Z' --output {output.zip} --silent"
        )
        output_folder = Path(output["shipdensity"]).parent
        shell("unzip {output.zip} -d {output_folder}")


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
    log:
        python="logs/build_region_shape/{region}.log",
        notebook="logs/build_region_shape/{region}.py.ipynb",
    notebook:
        "../actions/build_region_shape.py.ipynb"


rule build_cutout:
    message:
        "Downloading cutout for: {wildcards.region}."
    input:
        gpkg="resources/regions/{region}.gpkg",
    output:
        cutout="resources/cutouts/{region}.nc",
    params:
        era5_year=config["renewables"]["era5_year"],
    log:
        python="logs/build_cutout/{region}.log",
        notebook="logs/build_cutout/{region}.py.ipynb",
    notebook:
        "../actions/build_cutout.py.ipynb"


rule build_potentials_and_profiles:
    message:
        "Determining {wildcards.technology} potentials and generation profiles for {wildcards.region}."
    input:
        copernicus="resources/Copernicus_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif",
        gebco="resources/gebco/GEBCO_2021.nc",
        gebco_slope="resources/gebco/GEBCO_2021_slope.nc",
        wdpa="resources/WDPA_Feb2022.gpkg",
        wdpa_marine="resources/WDPA_WDOECM_Feb2022_marine.gpkg",
        shipping_routes="resources/shipdensity/shipdensity_global.tif",
        cutout="resources/cutouts/{region}.nc",
        region="resources/regions/{region}.gpkg",
    output:
        profiles="resources/profiles/{region}_{technology}.nc",
        area_mask="figures/{region}/{technology}/area_mask.png",
        capacity_factor_map="figures/{region}/{technology}/capacity_factors.html",
        potential_map="figures/{region}/{technology}/potential_map.html",
    params:
        technology_details=lambda w: config["renewables"][w.technology],
    wildcard_constraints:
        technology="(pvplant|wind_onshore|wind_offshore|csp_tower)",
    threads: 4
    log:
        python="logs/build_potentials_and_profiles/{region}_{technology}.log",
        notebook="logs/build_potentials_and_profiles/{region}_{technology}.py.ipynb",
    benchmark:
        "benchmarks/build_potentials_and_profiles/{region}_{technology}.csv"
    notebook:
        "../actions/build_potentials_and_profiles.py.ipynb"


# Convert atlite RES supply files into a single file
rule combine_atlite_supply:
    input:
        profiles=expand(
            "resources/profiles/{region}_{technology}.nc",
            technology=["wind_offshore", "wind_onshore", "pvplant"], # also supported: "csp-tower"
            allow_missing=True,
        ),
    output:
        supply="resources/supply_{region}.nc",
    log:
        python="logs/combine_atlite_supply/{region}.log",
        notebook="logs/combine_atlite_supply/{region}.py.ipynb",
    notebook:
        "../actions/combine_atlite_supply.py.ipynb"


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

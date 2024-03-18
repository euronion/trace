# SPDX-FileCopyrightText: 2021 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later

# - GEGIS rules: Require Julia to be setup seperately, see Readme.md - #

# Configure GlobalEnergyGIS to output files saves output files in a preconfigured location.
if config["GlobalEnergyGIS"].get("init_gegis", False) is True:

    rule set_GEGIS_base_dir:
        output:
            directory(config["GlobalEnergyGIS"]["base_dir"]),
        threads: 1
        script:
            "../actions/set_GEGIS_base_dir.jl"

    # Download auxiliary datasets for GEGIS
    rule download_GEGIS_dataset:
        output:
            Path(config["GlobalEnergyGIS"]["base_dir"]) / "protected.jld",  # not full list, only dependencies for rules below (proxy all others)
            Path(config["GlobalEnergyGIS"]["base_dir"]) / "gadm.tif",
        script:
            "../actions/download_GEGIS_datasets.jl"

    # Download ERA5 data for wind/solar/synthetic demand for GEGIS
    rule download_GEGIS_era5:
        output:
            Path(config["GlobalEnergyGIS"]["base_dir"]) / "era5wind{year}.h5",
            Path(config["GlobalEnergyGIS"]["base_dir"]) / "era5solar{year}.h5",
            Path(config["GlobalEnergyGIS"]["base_dir"]) / "era5temp{year}.h5",
        script:
            "../actions/download_GEGIS_era5.jl"


# Create region for GlobalEnergyGIS containing


# one or more areas defined by GADM (Database of Global Administrative Areas)
rule create_region:
    input:
        config["GlobalEnergyGIS"]["base_dir"] + "gadm.tif",
    output:
        config["GlobalEnergyGIS"]["base_dir"] + "regions_{region}.jld",
    script:
        "../actions/create_region.jl"


# Generate synthetic demand using machine learning / GDP / Pop projections
# under SSP scenarios with GlobalEnergyGIS
rule create_synthetic_demand:
    input:
        config["GlobalEnergyGIS"]["base_dir"] + "regions_{region}.jld",
        config["GlobalEnergyGIS"]["base_dir"] + "era5temp{era_year}.h5",
    output:
        config["GlobalEnergyGIS"]["base_dir"]
        + "output/SyntheticDemand_{region}_"
        + config["GlobalEnergyGIS"]["synthetic_demand"]["ssp_scenario"]
        + "-"
        + str(config["GlobalEnergyGIS"]["synthetic_demand"]["ssp_year"])
        + "_{era_year}.jld",
    script:
        "../actions/create_synthetic_demand.jl"


# Generate renewables potentials and time-series with GlobalEnergyGIS
rule create_renewables:
    input:
        config["GlobalEnergyGIS"]["base_dir"] + "era5wind{era_year}.h5",
        config["GlobalEnergyGIS"]["base_dir"] + "era5solar{era_year}.h5",
        config["GlobalEnergyGIS"]["base_dir"] + "regions_{region}.jld",
    output:
        wind=(
            config["GlobalEnergyGIS"]["base_dir"]
            + "output/GISdata_wind{era_year}_{region}.mat"
        ),
        solar=(
            config["GlobalEnergyGIS"]["base_dir"]
            + "output/GISdata_solar{era_year}_{region}.mat"
        ),
    script:
        "../actions/create_renewables.jl"


# Convert GEGIS file structure to netcdf for further processing
# (Supply: Renewables)
rule combine_GEGIS_supply:
    input:
        wind=(
            config["GlobalEnergyGIS"]["base_dir"]
            + "output/GISdata_wind{year}_{region}.mat"
        ),
        solar=(
            config["GlobalEnergyGIS"]["base_dir"]
            + "output/GISdata_solar{year}_{region}.mat"
        ),
    output:
        "resources/supply_{region}_{year}.nc",
    log:
        "logs/GEGIS/combine_GEGIS_supply_{region}_{year}.log",
    script:
        "../actions/combine_GEGIS_supply.py"


# Convert GEGIS file structure to netcdf for further processing
# (synthetic demand)
rule combine_GEGIS_demand:
    input:
        config["GlobalEnergyGIS"]["base_dir"]
        + "output/SyntheticDemand_{region}_"
        + config["GlobalEnergyGIS"]["synthetic_demand"]["ssp_scenario"]
        + "-"
        + str(config["GlobalEnergyGIS"]["synthetic_demand"]["ssp_year"])
        + "_{era_year}.jld",
    output:
        netcdf="resources/demand_{region}_{era_year}.nc",
        csv="resources/demand_annual_{region}_{era_year}.csv",
    log:
        "logs/GEGIS/combine_GEGIS_demand/{region}_{era_year}.log",
    script:
        "../actions/combine_GEGIS_demand.py"

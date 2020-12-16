configfile: "config.yaml"

rule create_additional_components:
    output:
        additional_components="resources/additional_components.pkl"
    log:
        python="logs/create_additional_components.log",
        notebook="logs/notebooks_processed/create_additional_components.ipynb"
    notebook:
        "actions/create_additional_components.py.ipynb"

rule create_network:
    input:
        efficiencies="data/efficiencies.csv",
        costs=f"../technology-data/outputs/costs_{config['scenario']['year']}.csv",
        wacc="data/wacc.csv",
        distances="data/distances.csv",
        network="escs/{esc}",
        additional_components="resources/additional_components.pkl"
    output:
        network=directory("resources/networks/{esc}/{from}-{to}")
    log:
        python="logs/create_network_{esc}_{from}-{to}.log",
        notebook="logs/notebooks_processed/create_network/{esc}/{from}-{to}.ipynb"
    notebook:
        "actions/create_network.py.ipynb"


# Allow for custom overwrite of annual electricity demand for exporters
demand_d = {
    "gegis": "resources/demand_TRACES_2013.nc",
    "custom": "data/overwrite/demand.csv"
}
demand_i = demand_d[config["scenario"]["synthetic_demand"].lower()]

rule attach_supply:
    input:
        supply="resources/supply_TRACES_2013.nc",
        demand=demand_i,
        costs=f"../technology-data/outputs/costs_{config['scenario']['year']}.csv",
        wacc="data/wacc.csv",
        network="resources/networks/{esc}/{from}-{to}",
        additional_components="resources/additional_components.pkl"
    output:
        network=directory("resources/networks_supplied/{esc}/{from}-{to}")
    log:
        python="logs/attach_supply_{esc}_{from}-{to}.log",
        notebook="logs/notebooks_processed/attach_supply/{esc}/{from}-{to}.ipynb"
    notebook:
        "actions/attach_supply.py.ipynb"
        
rule solve_network:
    input:
        network="resources/networks_supplied/{esc}/{from}-{to}",
        additional_components="resources/additional_components.pkl"
    output:
        network="results/{esc}/{from}-{to}/network.nc"
    log:
        python="logs/solve_network_{esc}_{from}-{to}.log",
        notebook="logs/notebooks_processed/solve_network/{esc}/{from}-{to}.ipynb"
    notebook:
        "actions/solve_network.py.ipynb"

## - GEGIS rules: Require Julia to be setup seperately, see Readme.md - ##

if config["GlobalEnergyGIS"].get("init_gegis", False) is True:
# Configure GlobalEnergyGIS to output files saves output files in a preconfigured location.
    rule set_GEGIS_base_dir:
        output:
            directory(config['GlobalEnergyGIS']['base_dir'])
        threads: 1
        script:
            "actions/set_GEGIS_base_dir.jl"

# Download auxiliary datasets for GEGIS
    rule download_GEGIS_dataset:
        output:
            # not full list, only dependencies for rules below (proxy all others)
            config['GlobalEnergyGIS']['base_dir']+"protected.jld",
            config['GlobalEnergyGIS']['base_dir']+"gadm.tif"
        script:
            "actions/download_GEGIS_datasets.jl"

# Download ERA5 data for wind/solar/synthetic demand for GEGIS
    rule download_GEGIS_era5:
        output:
            config['GlobalEnergyGIS']['base_dir']+"era5wind{year}.h5",
            config['GlobalEnergyGIS']['base_dir']+"era5solar{year}.h5",
            config['GlobalEnergyGIS']['base_dir']+"era5temp{year}.h5"
        script:
            "actions/download_GEGIS_era5.jl"

# Create region for GlobalEnergyGIS containing
# one or more areas defined by GADM (Database of Global Administrative Areas)
rule create_region:
    input:
        config['GlobalEnergyGIS']['base_dir']+"gadm.tif"
    output:
        config['GlobalEnergyGIS']['base_dir']+"regions_{region}.jld"
    script:
        "actions/create_region.jl"
        
# Generate synthetic demand using machine learning / GDP / Pop projections
# under SSP scenarios with GlobalEnergyGIS
rule create_synthetic_demand:
    input:
        config['GlobalEnergyGIS']['base_dir']+"regions_{region}.jld",
        config['GlobalEnergyGIS']['base_dir']+"era5temp{era_year}.h5"
    output:
        config['GlobalEnergyGIS']['base_dir']+"output/SyntheticDemand_{region}_"+config['GlobalEnergyGIS']['synthetic_demand']['ssp_scenario']+"-"+str(config['GlobalEnergyGIS']['synthetic_demand']['ssp_year'])+"_{era_year}.jld"
    script:
        "actions/create_synthetic_demand.jl"

# Generate renewables potentials and time-series with GlobalEnergyGIS
rule create_renewables:
    input:
       config['GlobalEnergyGIS']['base_dir']+"era5wind{era_year}.h5",
       config['GlobalEnergyGIS']['base_dir']+"era5solar{era_year}.h5",
       config['GlobalEnergyGIS']['base_dir']+"regions_{region}.jld"
    output:
        wind=config['GlobalEnergyGIS']['base_dir']+"output/GISdata_wind{era_year}_{region}.mat",
        solar=config['GlobalEnergyGIS']['base_dir']+"output/GISdata_solar{era_year}_{region}.mat"
    script:
        "actions/create_renewables.jl"
        
# Convert GEGIS file structure to netcdf for further processing
# (Supply: Renewables)
rule combine_GEGIS_supply:
    input:
        wind=config['GlobalEnergyGIS']['base_dir']+"output/GISdata_wind{year}_{region}.mat",
        solar=config['GlobalEnergyGIS']['base_dir']+"output/GISdata_solar{year}_{region}.mat"
    output:
        "resources/supply_{region}_{year}.nc"
    log:
        "logs/combine_GEGIS_supply_{region}_{year}.log"
    script:
        "actions/combine_GEGIS_supply.py"

# Convert GEGIS file structure to netcdf for further processing
# (synthetic demand)
rule combine_GEGIS_demand:
    input:
        config['GlobalEnergyGIS']['base_dir']+"output/SyntheticDemand_{region}_"+config['GlobalEnergyGIS']['synthetic_demand']['ssp_scenario']+"-"+str(config['GlobalEnergyGIS']['synthetic_demand']['ssp_year'])+"_{era_year}.jld"
    output:
        netcdf="resources/demand_{region}_{era_year}.nc",
        csv="resources/demand_annual_{region}_{era_year}.csv"
    log:
        "logs/combine_GEGIS_demand_{region}_{era_year}.log"
    script:
        "actions/combine_GEGIS_demand.py"
        
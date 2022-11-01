<!--
SPDX-FileCopyrightText: 2020-2022 Johannes Hampp
SPDX-License-Identifier: CC-BY-4.0
-->

[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/euronion/trace/main.svg)](https://results.pre-commit.ci/latest/github/euronion/trace/main) 
[![REUSE status](https://api.reuse.software/badge/github.com/euronion/trace/)](https://api.reuse.software/info/github.com/euronion/trace/)
[![Snakemake](https://img.shields.io/badge/snakemake-≥6.2.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

# TRACE - Transporting Renewables as Chemical Energy

## Setup and installation

1. Downloading the `git` repository:
   
    Download the main repository to somewhere (e.g. into a `projects` folder):

    ```
        (projects) > git clone https://github.com/euronion/trace.git
    ```

    This repository contains the main workflow, the energy supply chain definitions
    and all the necessary scripts tying it together.

2. Setup your `python` environment:

    a. Install the `python` environment manager `conda`.
    You need a `conda` compatible package manager installed.
    Recommended is [mamba](https://mamba.readthedocs.io/en/latest/).
    Alternatively use [miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/).

    b. Setup your `python` environment:

    ```
        conda env create -f envs/environment.yaml

        # or if you have mamba installed
        mamba env create -f envs/environment.yaml
    ```

    This installs all the packages and `python` dependencies.

3. Setup your `julia` environment (we used `julia-1.4`)

    a. Download and install `julia` following the installation instruction on the [Julia website](https://julialang.org/downloads/).

    b. Open a command shell and navigate into the `trace` folder
    ```
        cd trace
    ```
    c. Start `julia`
    ```
        (trace): julia
    ```
    d. Open the package manager of `julia`
    ```
        julia> ]
    ```
    e. Register and instantiate the environment
    ```
        (@v1.4) pkg> activate ./envs
        (trace) pkg> instantiate
    ```
4. (Optionatl) Setup a Gurobi solver license.
   For details and how to obtain a free academic license (if you are eligible), see the [documentation](https://gurobi.com/).
   If you don't want to use Gurobi you need to setup an [alternative linear solver compatible with PyPSA](https://pypsa.readthedocs.io/en/latest/installation.html#getting-a-solver-for-optimisation) and change the source code to be comptaible with the solver.
## Files and folders

The repository contains the follwing folder and files:

```
(trace): tree -L 2
.
├── actions     # Scripts for building and solving the model
│   ├── attach_supply.py.ipynb
│   ├── backup_run.py
│   ├── combine_GEGIS_demand.py
│   └── ...
├── config      # Input configuration for model and scenarios
│   ├── config.cern_link.yaml
│   ├── config.default.yaml
│   └── config.initial_paper.yaml
├── data        # Input data for the model
│   ├── technology-data
│   ├── distances.csv
│   ├── efficiencies.csv
│   ├── import_profiles
│   ├── overwrite
│   ├── shipping.csv
│   └── wacc.csv
├── envs        # Programming environment definitions for installing/reproducing the model
│   ├── environment.yaml
│   ├── environment.yaml.used
│   ├── JuliaManifest.toml
│   └── JuliaProject.toml
├── escs        # ESC definitions which serve as input; will be modified to create ESC-exp-imp specific models
│   ├── hvdc-to-elec
│   ├── hvdc-to-h2
│   ├── pipeline-ch4
│   └── ...
├── LICENSES    # License files for all licenses used in this repository
│   ├── CC0-1.0.txt
│   ├── CC-BY-4.0.txt
│   └── GPL-3.0-or-later.txt
├── logs        # (created automatically) log files for each model step, scenario and ESC
│   └── ...
├── README.md   # You are reading this
├── resources   # (created automatically) intermediary files created during model creation
│   └── ...
├── results     # (created automatically) Result files for each scenario/year/ESC/exporter-importer;
|   |           # Contains the optimised PyPSA networks and a CSV file with prominent results
│   ├── default
│   ├── lowhomogeneous
│   ├── results.csv
│   └── ...
├── rules       # Snakemake rule directory containing the workflow definitino for the model
│   ├── esc_construction.smk
│   ├── gegis.smk
│   ├── results.smk
│   └── solving.smk
└── Snakefile   # Snakemake file to combine all workflow files into the complete workflow
```

## Scenarios

This model uses *scenarios* to define a number of parameters identical across multiple model runs.
Scenarios are defined in `config/config.default.yaml` and `config/config.initial_paper.yaml`.
Identically named scenarios from the second file take precedence over scenarios from the first file.

Two central scenarios are the scenario `default` and `lowhomogeneous`.
Both are similar except for `default` assuming 10% WACC (`wacc: "homogeneous"`) and
`lowhomogeneous` assuming 5% WACC `wacc: "lowhomogeneous"`).
The values for `wacc: <value>` are references to column entries in `data/wacc.csv`.

## Energy Supply Chains (ESCs)

ESCs are defined in the `esc` subfolder as PyPSA networks.
For each ESC, its specific subfolder is parsed and converted into a region-specific PyPSA network representing the ESC.
Following files corresponding to PyPSA components are used for described purposes:

* 'loads.csv' : for implementing demand
    + Demand is implemented with a constant load (hourly average = annual total / 8760)
    + A store component (see later) is connected with losless links and without capital costs of arbitrary capacity
      to remove any influence by the time-dependency of the demand.
      With this approach the annual delivered energy becomes relevant, rather than satisfying the hourly demand in time.
      By convention, this store located on the import side and is named "buffer", is assigned with capital_cost=0.1 (to avoid numeric shennanigans).
* 'buses.csv' : basic building blocks to which all else attaches
    + buses represent different states of energy carriers or chemicals, e.g. CH4 (g), CH4 (l), CH4 compressed or water.
    + units of the buses are used for consitency checks and scaling costs/efficiencies in the model
    + by convention were energy is relevant, the bus carriers a power unit (preferrably MW); otherwise the bus uses m^3 or t (NOT m^3/h or t/h)
    + additional buses and links may be created to match the units in the costs input
    (e.g. converting t to m^3 for methanol storage, rather than changing the unit in the cost database)
* 'stores.csv' : Storing energy carriers or chemicals
* 'links.csv' : Converting between energy carriers, chemicals
    + Conversion defined via the 'data/efficiencies.csv' file
    + names of the links and connected buses must match to the efficiencies, except for a trailing (exp) or (imp)
    + names of the links and units of their bus0 or bus1 must match the names in the cost database, except for trailing (exp) or (imp) and different si prefixes for units
    + most buses are unidirectional with bus0 as input and bus1 as output; bus2, bus3 can be additional inputs
    + some buses are bidirectional with bus0 and bus 1 as input and output; bus2, bus3 not implemented for these types of links
    + special column 'scale_costs_based_on': This is an extension on the csv file. It specifies to which bus/energy carrier the 
        costs for this link in the cost data refers to. Values of this column are either 'bus0','bus1' or 'n/a' if no cost entry for this
        link exists because it is irrelevant.
* 'ships.csv' : Own extension for shipping based on PyPSA components
    + Contains properties of different types of ships
    + Data there must match the ship types in the cost database
    + for more details, see the 'Special case: Shipping' below
* Naming conventions
    + Trailing "(exp)" or "(imp)" are ignored during building of the network in bus/link/store names)
    + Also ignored are leading descriptors in the same braket, e.g. "(charging, exp)" or "(What do you hear?, imp")
    + Use "(exp)" and "(imp)" in names to indicate were a component is located - on the exporter or importer side of the model
* Costs
    + data taken from customised `technology-data` repository (CAPEX, FOM, lifetime)
    + capital costs p.a. are calculated using EAC
    + WACC assumptions in `config.yaml`
    + cost data is loaded during network creation (`actions/create_network.py.ipynb`)
    + adjustments to cost data (where necessary) are made right after network creation (`actions/create_network.py.ipynb`). The following adjustments are made:
        - Battery inverter cost: halfed; cost is for bidirectional inverter, model deploys two inverters with same power rating (enforced by constrained), thus doubling the cost for inverters
        - transmission distances: between regions are considered based on `data/distances.csv` and link capital costs per distance are scaled with the distance.

### Special case: Shipping
        - Shipping is special as it is a non-continuous transportation method
        - Shipping connections are defined with the non-standard PyPSA `ships.csv` in the respecitve ESC
        - `ships.csv` is read by `create_network.py.ipynb` and then translated into comparable PyPSA compoenents
        - Internally a ship is represented by a link for loading, a link for unloading and a link for representing energy loss during the journey and a store (cargo hold).
            The energy loss during the journey is either the energy required for propulsion or boil-off, whatever is larger.
            The idea is that boil-off may be used for propulsion as well and propulsion energy is always provided by the ship's cargo.
        - A delivery schedule is exogeneously provided, letting the loading and unloading links work with fixed delay representing the loading/unloading process and
            the travel time inbetween.
        - Instead of a ship with a fixed maximum capacity we think of it more as a group of ships = convoy which all travel together without an upper limit on capacity.
        - Additionally a constraint added in `solve_network.py.ipynb` for the solver ensures that for each convoy the annual loaded and unloaded amounts are equal.
        - For each shipping connection, multiple convoys are added such that all possible combinations are accounted for.
            In these combinations, the convoys can not load or unload while another convoy is loading/unloading.
        - The convoys are scheduled such that a as-constant-as-possible supply may be build by the model, i.e. a convoy arrives always shortly after the previous
            convoy left.
        - Travel time may be artificially increased to achieve a more constant supply by ship and avoid longer periods were not shipment arrives
            at the end of the year. The increase should be negiglible and in the range of a few hours per trip and is also logged for each scenario.
        - Additionally convoys may arrive a few hours later to fill gaps in the scheduling (again to smoothen the supply)
        - The model is provided with the investment costs per amount of cargo deliverd and may freely choose the amount to transport during for each shippment in each shipping "lane"
        - Limitation of the exogeneous scheduling:
          At the beginning of each scenario year there might be a larger supply gap of a few hours/days, where no ship arrives due to unfavourable combinations of
            journey time and loading/unloading duration
            
### Special case: LOHC

For LOHC transport the chemical costs also have to be taken into account.
This is hard-coded into the model in `actions/create_network.py.ipynb` in the method `override_costs_for_special_cases(n)`.
The following assumptions are made:
* DBT is used as LOHC
* Chemical costs are stored in cost database as 'LOHC chemical' in EUR/t
* The LOHC / H2 ratio is taken from efficiencies.csv , "LOHC hydrogenation"

The hard-coding modifies the following costs:
* The generators for LOHC chemicals: capital_costs become marginal_costs (per t of LOHC additionally required in the model to compensate losses)
* LOHC shipping lanes have an additional store for keeping track and transporting used, unloaded LOHC chemical

## Exporters and Importers

Exporting and importing regions are defined in `config/config.default.yaml` under the key
```
regions:
  TRACES:
    "AU": "Australia"
    "AR": "Argentina"
    "DE": "Germany"
    "DK": "Denmark"
    "EG": "Egypt"
    "ES": "Spain"
    "MA": "Morocco"
    "SA": "Saudi Arabia"
```

The model is setup to support imports from any of these regions to Germany (DE).
For each region a complete row entry in `data/wacc.csv` is necessary.
Further for each export-import relationship which is desired to be investigated, the different types of 
distances (as-the-crow-flies, as-the-hake-swims, sea route) are provided via `data/distances.csv`.

If a new region is supposed to be supported as importer or exporter, the new region has to be
included in all three locations with a complete entry or else the model will fail with an `KeyError`
when trying to run the model on a new region with incomplete entries.

## Running the model

> Note: This repository relies on `snakemake` as a workflow manager.
> For more information and usage of `snakemake` see the [documentation](https://snakemake.readthedocs.io/).

The model allows to simulate an energy supply chain (ESC) under various conditions.

## Technology data (costs)

Data technology costs (CAPEX, FOM, lifetime) are taken from the [technology-data GitHub repository](https://github.com/PyPSA/technology-data/)
and downloaded automatically into `data/technology-data/outputs/costs_<year>.csv`.
Depending on the `year` choosen for a scenario run (see below), the corresponding `costs_<year>.csv` file is used.
## Technology data (efficiencies)

Efficiencies for conversion processes, e.g. synthesis or hydrogen electrolysis are provided via `data/efficiencies.csv`.

### For a single case

A single model run or scenario run is a specific choice of the following:

* `scenario` (from the `config.default.yaml` and `config.initial_paper.yaml` files)
* energy supply chain (ESC) `esc`
* exporting region `exporter`
* importing region `importer`
* year (technology assumption) `year`

Use the following command where values `<...>` are substituted with the specific values:

```
    snakemake -jall results/<scenario>/<year>/<esc>/<exporter>-<importer>/results.csv
    
    # example
    snakemake -jall results/default/2030/shipping-lch4/AR-DE/results.csv
```

The command executes the necessary steps, solves the network(s) for all combinations of ESCs, exporting and importing countries,
extracts the results and merges the main results into a single `results/.../results.csv` file for manual or automatic result analysis.
In the same folder as the `results.csv` file there will also be saved a `network.nc` file.
This file is the solved, i.e. optimised, PyPSA network for the specific case modelled and can be opened using PyPSA

```
    import pypsa
    network = pypsa.Network(<path to network.nc file>)
```

All entries in `results.csv` are based on the `network.nc` and extracted for convinience using the Python notebook `actions/extract_results.py.ipynb`.

### Running multiple scenario runs

Using the `snakemake` functionality of parameter spaces one or scenarios runs can be defined and solved easily.
Each case to model and solve is specified by a row entry in `scenarios/default.csv`.
Every row must contain a valid entry for all 5 columns, else the workflow will fail.

Then all the scenario runs can be modelled and solved by calling

```
    snakemake -jall results/results.csv
```

This will create all the individual `results/.../results.csv` and `results/.../network.nc` files for each scenario run.
In addition it will create a file `results/results.csv` containing all individual results files combined into a single file for convinience.

## Sensitivity analysis

Sensitivity analysis is implemented through different scenarios, see `config/config.initial_paper.yaml` example
sensitivity scenarios and the `default` scenario in `config/config.default.yaml` for all implemented `modifiers` in the model.
Scenario runs for the sensitivity analysis are run as standard scenario cases, see the instructions above.

### Troubleshooting FAQ

* Workflow stops with "Address is already in use" error messages:
  In some cases scripts in the snakemake workflow can fail if too many `snakemake` rules are started simultaneously.
  This issue is related to the use of Jupyter notebook by the model.
  This problem cannot easily be avoided. Since the problem is not with the model but with the execution, an effective workaround
  is to tell `snakemake` to simply try to rerun the failed affected rules by adding  the `--restart-times <N>` flag added to call, e.g.

  ```
    snakemake -jall --restart-times 3 results/results.csv
  ```


* "Python kernel died" during `solve_network`:
   Memory consumption during solving step can exceeds the memory `snakemake` has allocated to the solving process.
   The model is setup to increase the reserved memory per repetition of the same `solve_network` step if the previous
   iteration has failed. Apply the same solution flags to `snakemake` as for the "Adress is already in use" problem bove.
   
  ```

## Conventions

* Chemicals are specified with their lower heating value (LHV) and use of LHV at this point is indicated
* Other input chemicals like water or CO2 are specified by their weight (t)
* bus units always represent a rate, i.e. "t/h" or "MW"
* The demand of the importing region is specified in `<esc>/loads.csv` in MWh/h, conversion table for hydrogen:

| energy carrier | hourly amount [t] | annual amount [t] | hourly amount [MWh] |annual amount [TWh] | HV assumed |
|----------------|-------------------|-------------------|---------------------|--------------------|------------|
| H2             | 410.96            | 3.6e6             | 13 698.65           | 120                | LHV        |

* "bus0" for links is always the input for which the cost.csv data is defined (preferred: electricity)
* "bus1" for links is always the output of the link
* RES are always attached to the bus "electricity (exp)" (in attach_supply.py.ipynb). This bus name must therefore exist in each ESC.
* In leap years (e.g. 2040) the leap day is ignored by dropping the 29th of February when encountered (= 365 days for all years)

## Modelling of renewable electricity supply and supply curves

* Supply time-series are synthetically generated (hourly curves, one representative year) using GlobalEnergyGIS (Julia package)
* Config for creating supply via `config/config.default.yaml`
* Supply is attached to networks via `actions/attach_supply.py.ipynb`
    1. LCoE for all different supply quality classes is determined
    2. Domestic demand is compared to LCoEs
    3. RES with lowest LCoE are reserved for domestic demand, such that the annual production from the RES classes matches the domestic demand.
    (no additional electricity demand or storage losses are assumed for satisfying the domestic electricity demand and no hourly profile,
    only the annual demand/supply).
    4. Residual RES are then attached to the network as PyPSA generators.
    
## Citing

This software is archived in Zenodo.
A permanent record of the most recent release version can be found here: [![DOI](https://zenodo.org/badge/317544698.svg)](https://zenodo.org/badge/latestdoi/317544698) .

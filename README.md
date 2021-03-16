<!--
SPDX-FileCopyrightText: 2020-2021 Johannes Hampp
SPDX-License-Identifier: CC-BY-4.0
-->

# TRACE - Transporting Renewables as Chemical Energy

## Setup and installation

1. Downloading `git` repositories:
   
    a. Download the main repository to somewhere (e.g. into a `projects` folder):

    ```
        (projects) > git clone https://github.com/euronion/trace.git
    ```

    This repository contains the main workflow, the energy supply chain definitions
    and all the necessary scripts tying it together.

    b. Download the `technology-data` repository into a parallel directory to the previous repository:

    ```
        (projects) > git clone -b new-techs/chemical-energy-carriers https://github.com/euronion/technology-data.git
    ```

    This repository contains the technology cost data.
    > Note:
    > That the repository branch is not `main` and the repository is currently **NOT** identical with https://github.com/PyPSA/technology-data/ .

    The resulting resulting structure should look like this

    ```
        (projects): tree -d -L 1 ./
        ./projects
        ├── technology-data
        ├── trace
        └── ...
    ```


2. Setup your `python` environment:

    a. Install the `python` environment manager `conda`.
    You need to either have [miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/) installed on your system.

    b. Setup your `python` environment:

    ```
        conda env create -f envs/environment.yaml
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

## Defining ESCs

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
            
### LOHC

For LOHC transport the chemical costs also have to be taken into account.
This is hard-coded into the model in `actions/create_network.py.ipynb` in the method `override_costs_for_special_cases(n)`.
The following assumptions are made:
* DBT is used as LOHC
* Chemical costs are stored in cost database as 'LOHC chemical' in EUR/t
* The LOHC / H2 ratio is taken from efficiencies.csv , "LOHC hydrogenation"

The hard-coding modifies the following costs:
* The generators for LOHC chemicals: capital_costs become marginal_costs (per t of LOHC additionally required in the model to compensate losses)
* LOHC shipping lanes have an additional store for keeping track and transporting used, unloaded LOHC chemical

### Conventions

* Chemicals are specified with their lower heating value (LHV) and use of LHV at this point is indicated
* Other input chemicals like water or CO2 are specified by their weight (t)
* The demand of the importing region is specified in `<esc>/loads.csv` in MWh/h

| energy carrier | hourly amount [t] | annual amount [t] | hourly amount [MWh] |annual amount [TWh] | HV assumed |
|----------------|-------------------|-------------------|---------------------|--------------------|------------|
| H2             | 410.96            | 3.6e6             | 13 698.65           | 120                | LHV        |

* "bus0" for links is always the input for which the cost.csv data is defined (preferred: electricity)
* "bus1" for links is always the output of the link
* RES are always attached to the bus "electricity (exp)" (in attach_supply.py.ipynb). This bus name must therefore exist in each ESC.
* In leap years (e.g. 2040) the leap day is ignored by dropping the 29th of February when encountered (= 365 days for all years)

## Supply

* Supply is generated (hourly curves, one year) using Global Energy GIS
* Config for creating supply via `config.yaml`
* Supply is attached to networks via `actions/attach_supply.py.ipynb`
    1. LCoE for all different supply quality classes is determined
    2. Domestic demand is compared to LCoEs
    3. RES with lowest LCoE are reserved for domestic demand, such that the annual production from the RES classes matches the domestic demand.
    (no additional electricity demand or storage losses are assumed for satisfying the domestic electricity demand and no hourly profile,
    only the annual demand/supply).
    4. Residual RES are then attached to the network as PyPSA generators.
    
## Solving

By default the commercial solver Gurobi is used.
Free academic licenses are available.
We use the 'barrier' solving algorithm and skip the crossover, 
as the additional accuracy from crossover is negligible for our problem.
We also make use of the "PreDual = 2" option from Gurobi,
which showed to dramatically increase solver speed for all ESCs involving long distance shipping.

As a fallback method the model automatically switches back to the default gurobi solver methods and
retries solving the scenario again.
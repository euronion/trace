# TRACE - Transporting Renewables as Chemical Energy

## Setup and installation

0. Download this repository to somewhere (into your current folder):

```
    git clone https://github.com/euronion/trace.git
```

1. Use `conda` for setting up environments.
2. Install `mamba` into a dedicated env.

```
    conda install -c conda-forge mamba
```

3. Install `snakemake` into a dedicated env.

```
    mamba create -c conda-forge -c bioconda -n snakemake snakemake
```

4. Setup `trace` environment from `environment.yaml` file

```
    conda env create -f environment.yaml
```
5. Download `technology-data` repository into a parallel folder to this repository.

```
    git clone https://github.com/PyPSA/technology-data.git
```
so that the resulting structure looks like this

```
    > tree -d -L 1
    .
    ├── technology-data
    ├── trace
    └── ...
```

## Defining ESCs

ESCs are defined in the `esc` subfolder as PyPSA networks.
For each ESC, its specific subfolder is parsed and converted into a region-specific PyPSA network representing the ESC.

* Implementation of demand
    + Demand is implemented with a constant load (hourly average = annual total / 8760)
    + A store is connected with losless links and without capital costs of arbitrary capacity
      to remove any influence by the time-dependency of the demand.
      With this approach the annual delivered energy becomes relevant, rather than satisfying the
      hourly demand in time.
* Costs
    + data taken from customised `technology-data` repository (CAPEX, FOM, lifetime)
    + capital costs p.a. are calculated using EAC
    + WACC assumptions in `config.yaml`
    + cost data is loaded during network creation (`actions/create_network.py.ipynb`)
    + adjustments to cost data (where necessary) are made right after network creation (`actions/create_network.py.ipynb`). The following adjustments are made:
        - Battery inverter cost: halfed; cost is for bidirectional inverter, model deploys two inverters with same power rating (enforced by constrained), thus doubling the cost for inverters
        - transmission distances: between regions are considered based on `data/distances.csv` and link capital costs per distance are scaled with the distance.
* Conversions
    + Conversions (e.g. electricity to hydrogen, seawater desalination) are all defined by PyPSA links via the `links.csv` in each respective ESC
* Special case of **shipping**
        - Shipping is special as it is a non-continuous transportation method
        - Shipping connections are defined with the non-standard PyPSA `ships.csv` in the respecitve ESC
        - `ships.csv` is read by `create_network.py.ipynb` and then translated into comparable PyPSA compoenents
        - Internally a ship is represented by a generator with negative generation in the exporting country and a generator with the same generation capacity but reduced by the
            shipping efficiency in its destination.
            A delivery schedule is exogeneously provided, letting the generators work with a fixed delay in between for the amount of time the loading or unloading process take.
        - Instead of a ship with a fixed maximum capacity we think of it more as a group of ships = convoy which all travel together without an upper limit on capacity.
        - Additionally a constraint added in `solve_network.py.ipynb` for the solver ensures that for each convoy the annual loaded and unloaded amounts are equal.
        - For each shipping connection, multiple convoys are added such that all possible combinations are accounted for.
            In these combinations, the convoys can not load or unload while another convoy is loading/unloading.
        - The model is provided with the investment costs per amount of cargo deliverd and may freely choose the amount to transport during for each shippment in each shipping "lane"
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
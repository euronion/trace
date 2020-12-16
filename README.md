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
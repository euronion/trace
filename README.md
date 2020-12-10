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
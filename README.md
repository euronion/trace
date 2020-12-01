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

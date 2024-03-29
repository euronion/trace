<!--
SPDX-FileCopyrightText: 2021 Johannes Hampp
SPDX-License-Identifier: CC-BY-4.0
-->

# Release Notes

## Upcoming release (TBD)

* New configuration option `import_profile`: The import profile for each ESC can now be set based on an
    external spreadsheet to distinguish between baseload or other types of import demand.
* New configuration option `import_buffer`: The import profile can now be buffered such that the import
    demand would need to be met e.g. as baseload (`import_buffer: False`) or on an annual basis 
    (`import_buffer: annually`).
    Buffering is disabled by default.
* Changed ESCs:
    + All explicit import buffering components have been removed from the ESCs (`stores.csv`) in favour 
      of the new configuration entry `import_buffer
* New ESCs, which end in electricity delivered rather than hydrogen (using an CCGT turbine at the ESC end):
    + HVDC-to-elec
    + pipeline-h2-to-elec
    + shipping-h2-to-elec
* Renames ESCs:
    + hvdc-to-h2 (formerly: "hvdc"): to reflect the delivered energy carrier
    (no change to the ESC itself)
* New WACC scenarios: 7.5% for all countries ("7.5% homogeneous")
* New configuration system:
    + config files are located in `config/` directory
    + `config.default.yaml` contains default configuration for the model
    + an additional config is used to overwrite default values and add scenario specific configuration.
      Set the appropriate `configfile:` path in `Snakefile.
* Include new import-exporter relationships:
    + Imports to CERN/Geneva (for CERN link idea)
    + Imports to various points in Europe from outside Europe (for PyPSA-EUR-Sec import exploration)
* Add scenario management via `scenarios/default.csv` instead of specifying combinations of scenario/ESC/year/exporter/importer via `config.yaml`s
* Add option for year-specific efficiencies: Efficiencies can now be specified in `data/efficiencies.csv` for either "all" years or year-specific values be given
  * Year-specific efficiencies for water electrolysis given, which affects most ESCs and thus was determined to be very important.
    All other efficiencies are constant over the years.

### Other:

* `pre-commit` continous integration used
* black coding style enforced
* Snakemake related files are formatted using `snakefmt`
* Fix: compatability of jupyter notebooks with new `nbconvert` versions
* `technology-data` (https://github.com/PyPSA/technology-data) is now automatically downloaded and included
  as `data/technology-data/<...>` files. The version being used by default is always the most up-to-date version
  and can be configured via `config/config.default.yaml` (new entry: `config["technology_data"]).

## Initial release (2021-08-03)

Initial release of the package, corresponding to the first of of [our preprint](https://arxiv.org/abs/2107.01092).

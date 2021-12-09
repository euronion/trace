<!--
SPDX-FileCopyrightText: 2021 Johannes Hampp
SPDX-License-Identifier: CC-BY-4.0
-->

# Release Notes

## Upcoming release (TBD)

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

### Other:

* `pre-commit` continous integration used
* black coding style enforced
* Snakemake related files are formatted using `snakefmt`

## Initial release (2021-08-03)

Initial release of the package, corresponding to the first of of [our preprint](https://arxiv.org/abs/2107.01092).
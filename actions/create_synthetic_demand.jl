# SPDX-FileCopyrightText: 2020-2021 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later

# Activate our julia environment
using Pkg
Pkg.activate("./envs")

using GlobalEnergyGIS

region = snakemake.wildcards["region"]
year = parse(Int64, snakemake.wildcards["era_year"])

ssp_scenario = snakemake.config["GlobalEnergyGIS"]["synthetic_demand"]["ssp_scenario"]
ssp_year = snakemake.config["GlobalEnergyGIS"]["synthetic_demand"]["ssp_year"]

predictdemand(gisregion=region, sspscenario=ssp_scenario, sspyear=ssp_year, era_year=year)

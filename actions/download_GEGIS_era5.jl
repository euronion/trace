# SPDX-FileCopyrightText: 2020-2021 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later

# Activate our julia environment
using Pkg
Pkg.activate("./envs")

using GlobalEnergyGIS


year = snakemake.wildcards["year"]


println("Downloading ERA5 data.")
era5download(year)

println("Converting and recompressing data.")
makewindera5(year=year)
makesolarera5(year=year)
maketempera5(year=year)

# Delete original files only after conversion/recompression was successfull
println("Deleting original ERA5 download files.")
clearvars_era5("wind", year=year)
clearvars_era5("solar", year=year)
clearvars_era5("temp", year=year)

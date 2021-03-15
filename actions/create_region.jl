# SPDX-FileCopyrightText: 2020-2021 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later

using GlobalEnergyGIS
using OrderedCollections

name = snakemake.wildcards["region"]
gadm_members = snakemake.config["regions"][name]
gadm_members = OrderedDict(gadm_members)
sort!(gadm_members)

# Convert to appropriate shape and call GADM function from GEGIS on regions
members = [[code GADM(area)] for (code, area) in gadm_members]
members = vcat(members...)

println(members)

saveregions(name, members)
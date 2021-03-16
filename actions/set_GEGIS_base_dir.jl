# SPDX-FileCopyrightText: 2020-2021 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later

#=
The handling of output locations in GEGIS is quite messy.
One defines the base directory for GEGIS into which
files are dumped based on a predefined folder and naming
scheme.
Setting the base directory location can only be done with
the command below which throws an error if the .cdsapi file
for configuration of the CDSAPI already exists, e.g. because
of another programm or GEGIS configuring it before.

Nevertheless subsequent calls to the function correctly set
the base dir correctly before the error is thrown.
I.e. "normal" behaviour is to see an error which we catch.
Only reaching beyong the function call constitutes an error
which we have to throw manually.
=#

# Activate our julia environment
using Pkg
Pkg.activate("./envs")


d = string(pwd(),"/",snakemake.config["GlobalEnergyGIS"]["base_dir"])

println("Creating new directory $d.")
mkpath(d)

try
    using GlobalEnergyGIS
    # Dummy data for CDSAPI - assume .cdsapi file already correctly create (no overwrite)
    saveconfig(d, 1, "1", agree_terms=true)
    
    # This consitutes an real error
    @assert(false, "This should never happen and indicates something went wrong.") 
catch
    # Expected (normal) behaviour
    println("GEGIS base directory set to: $d")
finally
    println("You can check ~/.GlobalEnergyGIS_config for correct settings.")
end

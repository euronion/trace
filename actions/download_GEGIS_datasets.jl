# SPDX-FileCopyrightText: 2020-2021 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later

using GlobalEnergyGIS

println("Downloading aux. datasets.")
download_datasets()

println("Rasterising datasets.")
rasterize_datasets(cleanup=:all)

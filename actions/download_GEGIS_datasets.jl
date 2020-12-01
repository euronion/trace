using GlobalEnergyGIS

println("Downloading aux. datasets.")
download_datasets()

println("Rasterising datasets.")
rasterize_datasets(cleanup=:all)

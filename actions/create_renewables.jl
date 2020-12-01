using GlobalEnergyGIS

region = snakemake.wildcards["region"]
year = parse(Int64, snakemake.wildcards["era_year"])
classes = snakemake.wildcards["classes"]

classes_solar_min = snakemake.config["GlobalEnergyGIS"]["solar"]["classes"][classes]["min"]
classes_solar_max = snakemake.config["GlobalEnergyGIS"]["solar"]["classes"][classes]["max"]

solar_options = GlobalEnergyGIS.solaroptions()
merge!(solar_options, Dict(Symbol(k) => v for (k,v) in snakemake.config["GlobalEnergyGIS"]["solar"]["options"]))
merge!(solar_options,
    Dict(:gisregion => region,
         :era_year => year,
         :pvclasses_min => classes_solar_min,
         :pvclasses_max => classes_solar_max,
         :cspclasses_min => classes_solar_min, 
         :cspclasses_max => classes_solar_max))

GISsolar(;solar_options...)

classes_wind_min = snakemake.config["GlobalEnergyGIS"]["wind"]["classes"][classes]["min"]
classes_wind_max = snakemake.config["GlobalEnergyGIS"]["wind"]["classes"][classes]["max"]

wind_options = GlobalEnergyGIS.windoptions()
merge!(wind_options, Dict(Symbol(k) => v for (k,v) in snakemake.config["GlobalEnergyGIS"]["wind"]["options"]))
merge!(wind_options,
    Dict(:gisregion => region,
         :era_year => year,
         :onshoreclasses_min => classes_wind_min,
         :onshoreclasses_max => classes_wind_max,
         :offshoreclasses_min => classes_wind_min,
         :offshoreclasses_max => classes_wind_max
        ))

GISwind(;wind_options...)
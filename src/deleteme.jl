cd("src")

using Random

include("../src/MimiMETA.jl")
include("../src/montecarlo.jl")

tp_used = "all"

Random.seed!(1)
model = full_model(;
        rcp="RCP4.5",
        ssp="SSP2",
        tdamage="pointestimate",
        slrdamage="mode",
        saf=ifelse(tp_used == "all" || occursin("SAF", tp_used), "Distribution mean", false),
        interaction=tp_used != "none",
        pcf=ifelse(tp_used == "all" || occursin("PCF", tp_used), "Fit of Hope and Schaefer (2016)", false),
        omh=ifelse(tp_used == "all" || occursin("OMH", tp_used), "Whiteman et al. beta 20 years", false),
        amaz=ifelse(tp_used == "all" || occursin("AMAZ", tp_used), "Cai et al. central value", false),
        gis=ifelse(tp_used == "all" || occursin("GIS", tp_used), "Nordhaus central value", false),
        wais=ifelse(tp_used == "all" || occursin("WAIS", tp_used), "Value", false),
        ism=ifelse(tp_used == "all" || occursin("ISM", tp_used), "Value", false),
        amoc=ifelse(tp_used == "all" || occursin("AMOC", tp_used), "IPSL", false))

myupdate_param!(model, :Consumption, :damagepersist, 1.)

run(model)

model[:Consumption, :T_AT][90]
model[:Consumption, :tempdamage][90, 185]
model[:Consumption, :tempdamage_bhm][90, 185]
model[:Consumption, :tempdamage_waid][90, 185]
model[:Consumption, :tempdamage_coacch][90, 185]

model[:Consumption, :T_country][90, 185]
model[:Consumption, :T_country_1990][185]



model[:Consumption, :tempdamage][50, 192]
model[:Consumption, :tempdamage_bhm][50, 192]
model[:Consumption, :tempdamage_waid][50, 192]
model[:Consumption, :T_AT][50]

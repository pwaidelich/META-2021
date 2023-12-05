using Random

# set working directory (adjust this to your system if required)
# cd("src")
include("../src/MimiMETA.jl")
include("../src/montecarlo.jl")


### create different parameter draws as returned by getsim()

# calculate MC draws incl. default persistence distribution (setting persist switch to true)
Random.seed!(1)
draws_with_persist = getsim(50, "Fit of Hope and Schaefer (2016)", # PCF
                                "Cai et al. central value", # AMAZ
                                "Nordhaus central value", # GIS
                                "Distribution", # WAIS
                                "Distribution", # SAF
                                true, # persit
                                false, # emuc
                                false) # prtp
draws_with_persist

# inspect persistence values in the draws object
draws_with_persist[!, :Consumption_damagepersist]

# create another draw object where persistence draws are overwritten with 0.
draws_with_persist_overwrite0p0 = copy(draws_with_persist)
draws_with_persist_overwrite0p0[!, :Consumption_damagepersist] .= 0.
draws_with_persist_overwrite0p0[!, :Consumption_damagepersist]

# create another draw object where the persistence column is dropped
# NOTE: this mimicks the result of setting the persistence switch in getsim() to false but ensures identical draws for other parameters
draws_without_persist = copy(draws_with_persist)
select!(draws_without_persist, Not(:Consumption_damagepersist))


### inspect Utility_world_disc_utility of the first MC draw in 2100

# a) persistence draws from Monte Carlo distribution (which should overwrite model default)
Random.seed!(1)
model = full_model(rcp="RCP4.5", ssp="SSP2")
results_with_persist = runsim(model, draws_with_persist,
                 true, # ism_used
                 true, # omh_used
                 true, # amoc_used
                 "Cai et al. central value", # AMAZ
                 "Distribution") # WAIS
results_with_persist[1][:Utility_world_disc_utility][91]
# result: -4.1535189484213926e7

# b) no persistence draws in the draws object (= model should default to 0.5)
Random.seed!(1)
model = full_model(rcp="RCP4.5", ssp="SSP2")
results_without_persist = runsim(model, draws_without_persist,
                 true, # ism_used
                 true, # omh_used
                 true, # amoc_used
                 "Cai et al. central value", # AMAZ
                 "Distribution") # WAIS
results_without_persist[1][:Utility_world_disc_utility][91]
# result: -4.1535189484213926e7

# c) persistence draws in the draws object overwritten with 0.0
Random.seed!(1)
model = full_model(rcp="RCP4.5", ssp="SSP2")
results_with_persist_overwrite0p0 = runsim(model, draws_with_persist_overwrite0p0,
                 true, # ism_used
                 true, # omh_used
                 true, # amoc_used
                 "Cai et al. central value", # AMAZ
                 "Distribution") # WAIS
results_with_persist_overwrite0p0[1][:Utility_world_disc_utility][91]
# result: -4.1535189484213926e7

# d) contrast this with results when we set the persistence manually to 0.0 in the model itself
Random.seed!(1)
model = full_model(rcp="RCP4.5", ssp="SSP2")
myupdate_param!(model, :Consumption, :damagepersist, 0.)
results_with_manualpersist0p0 = runsim(model, draws_without_persist,
                 true, # ism_used
                 true, # omh_used
                 true, # amoc_used
                 "Cai et al. central value", # AMAZ
                 "Distribution") # WAIS
results_with_manualpersist0p0[1][:Utility_world_disc_utility][91]
# result: -7.325960968365538e7

# ensure that 'Consumption_damagepersist' is the string used by getsim()
sum(names(draws_with_persist) .== "Consumption_damagepersist")
# check another example parameter
sum(names(draws_with_persist) .== "CO2Model_a0")

# check if model that has not been run yet has this parameter name
# NOTE: this is the if-condition in runsim() to update parameters
model = full_model(rcp="RCP4.5", ssp="SSP2")
inst = Mimi.build(model)
has_parameter(model.md, Symbol("Consumption_damagepersist"))
# -> returns false

# compare to example parameter
has_parameter(model.md, Symbol("CO2Model_a0"))
# -> returns true

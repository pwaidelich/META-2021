cd("src")
include("../src/MimiMETA.jl")
include("../src/montecarlo.jl")

# load JLD2 package for save_object and load_object functions
using JLD2
using Random

# dimensions are taken from src/basemodel.jl
countries = CSV.read("../data/pattern-scaling.csv", DataFrame).Country
years = collect(2010:2200)
seed_used = 25081776

global mc_samplesize_global = 2000

Random.seed!(seed_used)
global mc_draws_used = getsim(mc_samplesize_global,
                # NOTE: we use Hope & Schaefer PCF for consistency with Dietz et al 2021
                # this means that there is effectively no MC draws for the PCF component
                # alternatively, set this to 'Kessler probabilistic'
                "Fit of Hope and Schaefer (2016)", # PCF
                "Distribution", # AMAZ
                "Distribution", # GIS
                "Distribution",  # WAIS
                "Distribution", # SAF
               true, # persit
               false, # emuc
               false) # prtp

function run_equityvaluation_mc(mc_samplesize = 500; rcp = "RCP4.5", ssp = "SSP2", persist::Union{Float64, AbstractString} = "Distribution", seed_used = 25081776,
                                tp_used = "all", tdamage = "bhm_pointestimate", omh = "none",
                                dir_output = "C:/Users/pwaidelich/Downloads/GitHub - Local/equity_var_climate/data/META MC results")
    # check the function inputs
    if (isa(persist, AbstractString) && persist != "Distribution")
        error("persist can only be set to Distribution if it is a string")
    end

    if !(tp_used in ["all", "none", "PCF", "OMH", "AMAZ", "GIS", "WAIS", "ISM", "AMOC"])
        error("tp_used must be all, none, or the identifier of an individual TP (all-caps)")
    end

    if omh == "default"
        omhversion = "Whiteman et al. beta 20 years"
    elseif omh == "updated"
        omhversion = "Whiteman et al. beta 20 years updated"
    elseif omh == "none"
        omhversion = false
    else
        error("omh input not supported")
    end

    ## copy-pasted from src/MimiMETA.jl to have overview of function arguments:
    # function base_model(; rcp="CP-Base", ssp="SSP2", co2="Expectation", ch4="default", warming="Best fit multi-model mean", tdamage="none", slrdamage="none")
    # function full_model(; rcp="RCP4.5", ssp="SSP2", co2="Expectation", ch4="default", warming="Best fit multi-model mean", tdamage="pointestimate", slrdamage="mode", saf="Distribution mean", interaction=true, pcf="Fit of Hope and Schaefer (2016)", omh="Whiteman et al. beta 20 years", amaz="Cai et al. central value", gis="Nordhaus central value", wais="Value", ism="Value", amoc="IPSL", nonmarketdamage=false)

    Random.seed!(seed_used)
    global model = full_model(;
        rcp=rcp,
        ssp=ssp,
        tdamage=tdamage,
        slrdamage="mode",
        # NOTE: we always use SAF
        saf="Distribution mean",
        interaction=tp_used != "none",
        pcf=ifelse(tp_used == "all" || occursin("PCF", tp_used), "Fit of Hope and Schaefer (2016)", false),
        omh=ifelse(tp_used == "all" || occursin("OMH", tp_used), omhversion, false),
        amaz=ifelse(tp_used == "all" || occursin("AMAZ", tp_used), "Cai et al. central value", false),
        gis=ifelse(tp_used == "all" || occursin("GIS", tp_used), "Nordhaus central value", false),
        wais=ifelse(tp_used == "all" || occursin("WAIS", tp_used), "Value", false),
        ism=ifelse(tp_used == "all" || occursin("ISM", tp_used), "Value", false),
        amoc=ifelse(tp_used == "all" || occursin("AMOC", tp_used), "IPSL", false))

    ## Update persistence parameter phi (hard-coded default: 0.5; MC default: uniform from 0 to 1)
    mc_draws_mutated = copy(mc_draws_used)

    if !isa(persist, AbstractString)
        mc_draws_mutated.Consumption_damagepersist .= persist
    end

    Random.seed!(seed_used)
    results = runsim(model, mc_draws_mutated,
                     tp_used == "all" || occursin("ISM", tp_used), # ism_used
                     omhversion != false && (tp_used == "all" || occursin("OMH", tp_used)), # omh_used
                     tp_used == "all" || occursin("AMOC", tp_used), # amoc_used
                     true, # saf_used (always true)
                     ifelse(tp_used == "all" || occursin("AMAZ", tp_used), "Distribution", "none"), # AMAZ
                     ifelse(tp_used == "all" || occursin("WAIS", tp_used), "Distribution", "none") # WAIS
                     ) # save all random variables)

    # save out as a pickle-like HDF object
    save_object(string(dir_output, "/jld2_files/", mc_samplesize,
                        # remove special characters from scenario names that can cause filename issues: ".", "/"
                        "_", replace(replace(rcp, "." => ""), "/" => ""), replace(ssp, "." => ""),
                        "_persist", persist,
                        "_tip", tp_used, "_tdamage", tdamage,
                        "_omh", omh,
                        ".jld2"),
                results)

    # return the results
    results
end

# test the function quickly
# temp = run_equityvaluation_mc(50, rcp = "RCP4.5", ssp = "SSP2", persist = "Distribution", tp_used = "all")

# procedure for nyears X ncountries matrices
function convert_mcresults_to_longdf(results, variable)
    if typeof(results[1][variable]) == Matrix{Union{Missing, Float64}} && size(results[1][variable]) == (length(years), length(countries))
        df_results = hcat(DataFrame(mc_run = 1, year = years), DataFrame(results[1][variable], countries))

        for jj in 2:length(results)
            df_results = vcat(df_results, hcat(DataFrame(mc_run = jj, year = years), DataFrame(results[jj][variable], countries)))
        end
    else
        # condense results from MC runs into a single DataFrame
        df_results = DataFrame(mc_run = Int64[], variable = Symbol[], value = Float64[])

        for jj in 1:length(results)
            push!(df_results, [jj, variable, results[jj][variable]])
        end
    end

    df_results
end

# test the function
# convert_mcresults_to_longdf(temp, :Consumption_lossfactor_temp)

function extract_lossfactors(mc_samplesize=nothing; rcp=nothing, ssp=nothing, persist=nothing, tp_used = "all",
                             tdamage = "bhm_pointestimate", omh = "none",
                             dir_output = "C:/Users/pwaidelich/Downloads/GitHub - Local/equity_var_climate/data/META MC results")

    # combine the parameter combinations into an identifier string
    parameter_id = string("n", mc_samplesize, "_", replace(replace(rcp, "." => ""), "/" => ""), replace(ssp, "." => ""),
                          "_persist", persist, "_tip", tp_used, "_tdamage", tdamage, "_omh", omh)

    # print out the parameters to track progress
    print(parameter_id)

    # if the file already exists, skip this interation
    if isfile(string(dir_output, "/lossfactors/lossfactor_conspc/lossfactor_conspc_", parameter_id, ".csv"))
        print("File already exists. Skipping...")
        return "Done"
    end

    # run equity valuations with the set parameters
    results = run_equityvaluation_mc(mc_samplesize,
                                     rcp =  rcp,
                                     ssp = ssp,
                                     tp_used = tp_used,
                                     persist = persist,
                                     tdamage = tdamage,
                                     omh = omh)

    # loop through different damage factors of interest
    for lossfactor_name_jj in ["lossfactor_conspc"] #, "lossfactor_persistence",
                                #"lossfactor_SLR", "lossfactor_extradamages", "lossfactor_temp"]
        # track Progress
        print(string("Extracting ", lossfactor_name_jj, "..."))

        # extract the damage factor of interest
        df_gdppc_damagefactor = convert_mcresults_to_longdf(results, Symbol(string("Consumption_", lossfactor_name_jj)))

        # add meta information
        insertcols!(df_gdppc_damagefactor, :rcp => rcp)
        insertcols!(df_gdppc_damagefactor, :ssp => ssp)
        insertcols!(df_gdppc_damagefactor, :tp_used => tp_used)
        insertcols!(df_gdppc_damagefactor, :persist => persist)
        insertcols!(df_gdppc_damagefactor, :tdamage => tdamage)

        # export as CSV
        CSV.write(string(dir_output, "/lossfactors/", lossfactor_name_jj, "/",
        lossfactor_name_jj, "_", parameter_id, ".csv"), df_gdppc_damagefactor)
    end

    return "Done"
end

# test the function
# extract_lossfactors(50, rcp = "RCP4.5", ssp = "SSP2", persist = "Distribution",
#                    tp_used = "all")

# run MC across scenarios & relevant settings, saving out results as .jld2 files
for omh_jj in ["none"]#, "default"]
    for persist_jj in [1.] #"Distribution"] #, 0., 1.] # 1.0 = no persistence, 0.0 = full persistence
        for tdamage_jj in ["waid_pointestimate", "coacch_central"] # "bhm_distribution"] #, "waid_distribution"] #, "coacch_central"]
            for (rcp_jj, ssp_jj) in [("RCP4.5", "SSP2")] #("RCP3-PD/2.6", "SSP2"), ("RCP4.5", "SSP2"), ("RCP8.5", "SSP2")]
                for tp_used_jj in ["all", "none", "PCF", "GIS", "WAIS", "AMAZ", "OMH", "ISM", "AMOC"]

                    # skip parameter combinations that are not of interest
                    #if (rcp_jj, ssp_jj) != ("RCP4.5", "SSP2") && (!(tp_used_jj in ["all", "none"]) || !(tdamage_jj in ["coacch_central", "waid_pointestimate"]) || persist_jj != 1. || omh_jj != "none")
                    #    continue
                    #elseif persist_jj != "1.0" && (!(tp_used_jj in ["all", "none"]) || !(tdamage_jj in ["bhm_pointestimate", "waid_pointestimate", "coacch_central"]) || omh_jj != "none")
                    #    continue
                    #elseif !(tp_used_jj in ["all", "none"]) && (tdamage_jj != "waid_pointestimate" || omh_jj != "none")
                    #    continue
                    #elseif tdamage_jj != "waid_pointestimate" && omh_jj != "none"
                    #    continue
                    #elseif omh_jj != "none" && tp_used_jj == "none"
                    #    continue
                    #end

                    extract_lossfactors(mc_samplesize_global, rcp=rcp_jj, ssp=ssp_jj, persist=persist_jj, tp_used=tp_used_jj,
                                        tdamage = tdamage_jj, omh=omh_jj)
                end
            end
        end
    end
end

# do some runs manually
extract_lossfactors(mc_samplesize_global, rcp="RCP8.5", ssp="SSP5", persist=1., tp_used="all",
                    tdamage = "coacch_central", omh="none")

#### cost to investors
# map the directory
dir_meta_results = "C:/Users/pwaidelich/Downloads/GitHub - Local/equity_var_climate/data/META MC results"

# baseline without pulse
Random.seed!(seed_used)
model = full_model(;
    rcp="RCP4.5",
    ssp="SSP2",
    tdamage="pointestimate",
    slrdamage="mode",
    saf="Distribution mean",
    interaction=true,
    pcf="Fit of Hope and Schaefer (2016)",
    omh="Whiteman et al. beta 20 years",
    amaz="Cai et al. central value",
    gis="Nordhaus central value",
    wais="Value",
    ism="Value",
    amoc="IPSL")

run(model)

df_lossfactor_base = hcat(DataFrame(year = years), DataFrame(model[:Consumption, :lossfactor_conspc], countries))
CSV.write(string(dir_meta_results, "/lossfactors/lossfactor_conspc_pulses/singlerun_baseline.csv"), df_lossfactor_base)

pulse_size_co2 = 1
# IPCC AR6 Table 7.15: 27 g CO2e GWP100 for non-fossil CH4, 29.8 for fossil CH4
# NOTE: ch4_extra is in Mt, not Gt
pulse_size_ch4 = 1*(1/27)*1000
pulse_year = 2024
pulse_index = findfirst(dim_keys(model, :time) .== pulse_year)

co2_extra_new = zeros(dim_count(model, :time))
co2_extra_new[pulse_index] = pulse_size_co2

ch4_extra_new = zeros(dim_count(model, :time))
ch4_extra_new[pulse_index] = pulse_size_ch4

# CO2 pulse
Random.seed!(seed_used)
model_co2pulse = full_model(;
    rcp="RCP4.5",
    ssp="SSP2",
    tdamage="pointestimate",
    slrdamage="mode",
    saf="Distribution mean",
    interaction=true,
    pcf="Fit of Hope and Schaefer (2016)",
    omh="Whiteman et al. beta 20 years",
    amaz="Cai et al. central value",
    gis="Nordhaus central value",
    wais="Value",
    ism="Value",
    amoc="IPSL")

myupdate_param!(model_co2pulse, :CO2Model, :co2_extra, co2_extra_new)
run(model_co2pulse)

df_lossfactor_co2pulse = hcat(DataFrame(year = years), DataFrame(model_co2pulse[:Consumption, :lossfactor_conspc], countries))
CSV.write(string(dir_meta_results, "/lossfactors/lossfactor_conspc_pulses/singlerun_co2pulse.csv"), df_lossfactor_co2pulse)

# CH4 pulse_size
Random.seed!(seed_used)
model_ch4pulse = full_model(;
    rcp="RCP4.5",
    ssp="SSP2",
    tdamage="pointestimate",
    slrdamage="mode",
    saf="Distribution mean",
    interaction=true,
    pcf="Fit of Hope and Schaefer (2016)",
    omh="Whiteman et al. beta 20 years",
    amaz="Cai et al. central value",
    gis="Nordhaus central value",
    wais="Value",
    ism="Value",
    amoc="IPSL")

myupdate_param!(model_ch4pulse, :CH4Model, :ch4_extra, ch4_extra_new)
run(model_ch4pulse)

df_lossfactor_ch4pulse = hcat(DataFrame(year = years), DataFrame(model_ch4pulse[:Consumption, :lossfactor_conspc], countries))
CSV.write(string(dir_meta_results, "/lossfactors/lossfactor_conspc/lossfactor_conspc_singlerun_ch4pulse_n1_RCP45SSP2_persistDistribution_tipall_tdamageBHM.csv"), df_lossfactor_ch4pulse)

# placebo
Random.seed!(seed_used)
model_placebo = full_model(;
    rcp="RCP4.5",
    ssp="SSP2",
    tdamage="pointestimate",
    slrdamage="mode",
    saf="Distribution mean",
    interaction=true,
    pcf="Fit of Hope and Schaefer (2016)",
    omh="Whiteman et al. beta 20 years",
    amaz="Cai et al. central value",
    gis="Nordhaus central value",
    wais="Value",
    ism="Value",
    amoc="IPSL")

myupdate_param!(model_placebo, :CH4Model, :ch4_extra, zeros(dim_count(model, :time)))
run(model_placebo)

df_lossfactor_placebo = hcat(DataFrame(year = years), DataFrame(model_placebo[:Consumption, :lossfactor_conspc], countries))
CSV.write(string(dir_meta_results, "/lossfactors/lossfactor_conspc_pulses/singlerun_placebo.csv"), df_lossfactor_placebo)

cd("src")
include("../src/MimiMETA.jl")
include("../src/montecarlo.jl")

# load JLD2 package for save_object and load_object functions
using JLD2
using Random

# dimensions are taken from src/basemodel.jl
countries = CSV.read("../data/pattern-scaling.csv", DataFrame).Country
years = collect(2010:2200)

function run_equityvaluation_mc(mc_samplesize = 50; rcp = "RCP4.5", ssp = "SSP2", persist::Union{Float64, AbstractString} = "Distribution", seed_used = 25081776,
                                tp_used = "all",
                                dir_output = "C:/Users/pwaidelich/Downloads/GitHub - Local/equity_var_climate/analysis")
    ##  Load model and select configuration
    Random.seed!(seed_used)

    # check the function inputs
    if (isa(persist, AbstractString) && persist != "Distribution")
        error("persist can only be set to Distribution if it is a string")
    end

    if !(tp_used in ["all", "none", "SAF", "PCF", "OMH", "AMAZ", "GIS", "WAIS", "ISM", "AMOC"])
        error("tp_used must be all, none, or the identifier of an individual TP (all-caps)")
    end

    ## copy-pasted from src/MimiMETA.jl to have overview of function arguments:
    # function base_model(; rcp="CP-Base", ssp="SSP2", co2="Expectation", ch4="default", warming="Best fit multi-model mean", tdamage="none", slrdamage="none")
    # function full_model(; rcp="RCP4.5", ssp="SSP2", co2="Expectation", ch4="default", warming="Best fit multi-model mean", tdamage="pointestimate", slrdamage="mode", saf="Distribution mean", interaction=true, pcf="Fit of Hope and Schaefer (2016)", omh="Whiteman et al. beta 20 years", amaz="Cai et al. central value", gis="Nordhaus central value", wais="Value", ism="Value", amoc="IPSL", nonmarketdamage=false)

    global model = full_model(;
        rcp=rcp,
        ssp=ssp,
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

    ## Update persistence parameter phi (hard-coded default: 0.5; MC default: uniform from 0 to 1)
    if !isa(persist, AbstractString)
        update_param!(model, :Consumption, :damagepersist, persist)
    end

    Random.seed!(seed_used)
    draws = getsim(mc_samplesize, ifelse(tp_used == "all" || occursin("PCF", tp_used), "Kessler probabilistic", "none"), # PCF
                   ifelse(tp_used == "all" || occursin("AMAZ", tp_used), "Distribution", "none"), # AMAZ
                   ifelse(tp_used == "all" || occursin("GIS", tp_used), "Distribution", "none"), # GIS
                   ifelse(tp_used == "all" || occursin("WAIS", tp_used), "Distribution", "none"), # WAIS
                   ifelse(tp_used == "all" || occursin("SAF", tp_used), "Distribution", "none"), # SAF
                   # if persist is 'Distribution' (and hence a string), we set the persit argument to true to draw from the MC distribution
                   isa(persist, AbstractString), # persit
                   false, # emuc
                   false) # prtp

    results = runsim(model, draws, tp_used == "all" || occursin("ISM", tp_used), # ism_used
                                   tp_used == "all" || occursin("OMH", tp_used), # omh_used
                                   tp_used == "all" || occursin("AMOC", tp_used), # amoc_used
                                   ifelse(tp_used == "all" || occursin("AMAZ", tp_used), "Distribution", "none"), # AMAZ
                                   ifelse(tp_used == "all" || occursin("WAIS", tp_used), "Distribution", "none") # WAIS
                                   ) # save all random variables)

    # save out as a pickle-like HDF object
    save_object(string(dir_output, "/lossfactors/jld2_files/", mc_samplesize,
                        # remove special characters from scenario names that can cause filename issues: ".", "/"
                        "_", replace(replace(rcp, "." => ""), "/" => ""), replace(ssp, "." => ""),
                        "_persist", persist,
                        "_tip", tp_used, ".jld2"),
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
                             dir_output = "C:/Users/pwaidelich/Downloads/GitHub - Local/equity_var_climate/analysis")

    # combine the parameter combinations into an identifier string
    parameter_id = string("n", mc_samplesize, "_", replace(replace(rcp, "." => ""), "/" => ""), replace(ssp, "." => ""),
                          "_persist", persist, "_tip", tp_used)

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
                                     persist = persist)

    # loop through different damage factors of interest
    for lossfactor_name_jj in ["lossfactor_conspc", "lossfactor_persistence",
                                "lossfactor_SLR", "lossfactor_extradamages", "lossfactor_temp"]
        # track Progress
        print(string("Extracting ", lossfactor_name_jj, "..."))

        # extract the damage factor of interest
        df_gdppc_damagefactor = convert_mcresults_to_longdf(results, Symbol(string("Consumption_", lossfactor_name_jj)))

        # add meta information
        insertcols!(df_gdppc_damagefactor, :rcp => rcp)
        insertcols!(df_gdppc_damagefactor, :ssp => ssp)
        insertcols!(df_gdppc_damagefactor, :tp_used => tp_used)
        insertcols!(df_gdppc_damagefactor, :persist => persist)

        # export as CSV
        CSV.write(string(dir_output, "/lossfactors/", lossfactor_name_jj, "/",
        lossfactor_name_jj, "_", parameter_id, ".csv"), df_gdppc_damagefactor)
    end
end

# test the function
# extract_lossfactors(50, rcp = "RCP4.5", ssp = "SSP2", persist = "Distribution",
#                    tp_used = "all")

# set the overall Monte Carlo sample size
mc_samplesize = 50

# run MC across scenarios & relevant settings, saving out results as .jld2 files
for (rcp_jj, ssp_jj) in [("RCP3-PD/2.6", "SSP2"), ("RCP4.5", "SSP2"), ("RCP8.5", "SSP2"),
                         ("RCP3-PD/2.6", "SSP5"), ("RCP4.5", "SSP5"), ("RCP8.5", "SSP5")]
    for persist_jj in ["Distribution", 0., 1., 0.5]
        for tp_used_jj in ["all", "none", "SAF", "PCF", "GIS", "WAIS", "AMAZ", "OMH", "ISM", "AMOC"]

            # skip parameter combinations that are not of interest
            if !(tp_used_jj in ["all", "none"]) && ((rcp_jj, ssp_jj) != ("RCP4.5", "SSP2") || persist_jj != "Distribution")
                continue
            elseif persist_jj != "Distribution" && ((rcp_jj, ssp_jj) != ("RCP4.5", "SSP2"))
                continue
            end

            extract_lossfactors(mc_samplesize, rcp=rcp_jj, ssp=ssp_jj, persist=persist_jj, tp_used=tp_used_jj)
        end
    end
end

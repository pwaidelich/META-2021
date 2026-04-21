# set the working directory
cd("src")

# execute the required META scripts
include("../src/MimiMETA.jl")
include("../src/montecarlo.jl")

# load JLD2 package for save_object and load_object functions
using JLD2
using Random
# load DataFrames and Threads for multi-threading
using DataFrames
using Base.Threads

# extract country ISO3 codes from META
countries = CSV.read("../data/pattern-scaling.csv", DataFrame).Country

# load the ISO3 codes of countries that feature in MSCI indices
countries_with_marketcap = CSV.read("../data/META_countries_with_marketcap.csv", DataFrame).iso3

# map the years in META manually
years = collect(2010:2200)

# fix the master seed, the MC sample size and the output directory globally
global seed_global = 25081776
global mc_samplesize_global = 10000
# global dir_output_global = "C:/Users/pwaidelich/Downloads/GitHub - Local/equity_var_climate/data/META MC results"
global dir_output_global = "Y:/gem/10 Research/10 Projects/GREENFIN/70 Equity VaR Tipping Points/META_results_Aug2024"
# global dir_output_global = "//d/groups/gess/cfp/gem/10 Research/10 Projects/GREENFIN/70 Equity VaR Tipping Points/META_results_Aug2024"

#####################################################
############## DRAW MONTE CARLO PARAMETERS ##########
#####################################################

# draw the MC parameters via getsim()
Random.seed!(seed_global)
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

### add the damage function draws
# NOTE: we do this outside of Mimi's functionality to facilitate multivariate distributions
# first, calibrate multivariate Gaussian using point estimates & vcov matrix from BHM, main spec
bhm_mvn = MvNormal([0.012718353, -0.00048709], [0.0000143458 -0.000000375824; -0.000000375824 0.0000000140157])
# then, draw parameters (jointly)
Random.seed!(seed_global)
bhm_draws = rand(bhm_mvn, mc_samplesize_global)
mc_draws_used.Consumption_beta1_bhm = bhm_draws[1, :]
mc_draws_used.Consumption_beta2_bhm = bhm_draws[2, :]

# set the coacch seeds to integers starting from 1 and ending at MC sample size
mc_draws_used.Consumption_seed_coacch = collect(1:mc_samplesize_global)

# export the MC draws used as a CSV file for reproducibility
CSV.write(string(dir_output_global, "/mc_draws_used_n", mc_samplesize_global, ".csv"), mc_draws_used)


#####################################################
############ DEFINE MC ANALYSIS FUNCTION ############
#####################################################

# write a function that runs META's MC mode with the required specification
function run_equityvaluation_mc(mc_samplesize = mc_samplesize_global; rcp = "RCP4.5", ssp = "SSP2", persist::Union{Float64, AbstractString} = 1., seed_used = seed_global,
                                tp_used = "all", tdamage = "coacch_distribution", omh = "none", amoc = "none", saf = "Distribution mean",
                                dir_output = dir_output_global)
    # check the function inputs
    if mc_samplesize > size(mc_draws_used)[1]
        error("Selected mc_samplesize exceeds the no. of rows in mc_draws_used")
    end

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

    if amoc == "none"
        amocversion = false
    elseif amoc == "Cai"
        amocversion = "Hadley"
    else
        amocversion = amoc
    end

    if !(amoc in ["IPSL", "BCM", "HADCM", "Hadley", "none", "Cai"])
        error("amoc input not supported")
    end

    # if saf is 'none', we set it to false
    if saf == "none"
        safversion = false
    else
        safversion = saf
    end

    if !(saf in ["Distribution mean", "none"])
        error("saf input not supported")
    end

    ## copy-pasted from src/MimiMETA.jl to have overview of function arguments:
    # function base_model(; rcp="CP-Base", ssp="SSP2", co2="Expectation", ch4="default", warming="Best fit multi-model mean", tdamage="none", slrdamage="none")
    # function full_model(; rcp="RCP4.5", ssp="SSP2", co2="Expectation", ch4="default", warming="Best fit multi-model mean", tdamage="pointestimate", slrdamage="mode", saf="Distribution mean", interaction=true, pcf="Fit of Hope and Schaefer (2016)", omh="Whiteman et al. beta 20 years", amaz="Cai et al. central value", gis="Nordhaus central value", wais="Value", ism="Value", amoc="IPSL", nonmarketdamage=false)

    # set up the model
    Random.seed!(seed_used)
    global model = full_model(;
        rcp=rcp,
        ssp=ssp,
        tdamage=tdamage,
        slrdamage="mode",
        saf=safversion, # this will be false if saf is 'none'
        interaction=tp_used != "none",
        pcf=ifelse(tp_used == "all" || occursin("PCF", tp_used), "Fit of Hope and Schaefer (2016)", false),
        omh=ifelse(tp_used == "all" || occursin("OMH", tp_used), omhversion, false),
        amaz=ifelse(tp_used == "all" || occursin("AMAZ", tp_used), "Cai et al. central value", false),
        gis=ifelse(tp_used == "all" || occursin("GIS", tp_used), "Nordhaus central value", false),
        wais=ifelse(tp_used == "all" || occursin("WAIS", tp_used), "Value", false),
        ism=ifelse(tp_used == "all" || occursin("ISM", tp_used), "Value", false),
        amoc=ifelse(tp_used == "all" || occursin("AMOC", tp_used), amocversion, false))

    # update persistence parameter phi (hard-coded default: 0.5; MC default: uniform from 0 to 1)
    mc_draws_mutated = copy(mc_draws_used)
    if !isa(persist, AbstractString)
        mc_draws_mutated.Consumption_damagepersist .= persist
    end

    # if we use the AMOC collapse a la Cai et al. (2016), overwrite the model parameter (defaults to 0.)
    if amoc == "Cai"
        myupdate_param!(model, :AMOC, :use_AMOC_Cai_impacts, 1.)
    end

    # create an identifier string for the model settings used
    parameter_id = string("n", mc_samplesize, "_", replace(replace(rcp, "." => ""), "/" => ""), replace(ssp, "." => ""),
                          "_persist", persist, "_tip", tp_used, "_tdamage", tdamage, "_omh", omh, "_amoc", amoc)

    # if SAF is deactivated, add this to the identifier string
    if saf == "none"
        parameter_id = string(parameter_id, "_safnone")
    end

    # print out for tracking progress
    print(parameter_id)

    # map (& if required, create) the subfolder where GDP loss CSV files will be stored
    subfolder_path = string(dir_output, "/", parameter_id)
    if !isdir(subfolder_path)
        mkdir(subfolder_path)
    end

    if !isdir(string(subfolder_path, "/lossfactor_conspc"))
        mkdir(string(subfolder_path, "/lossfactor_conspc"))
    end

    if !isdir(string(subfolder_path, "/T_AT"))
        mkdir(string(subfolder_path, "/T_AT"))
    end

    if !isdir(string(subfolder_path, "/SLR"))
        mkdir(string(subfolder_path, "/SLR"))
    end

    # if the N-th MC run in the subfolder already exists (= spec has already been run), exit
    if isfile(string(subfolder_path, "/lossfactor_conspc/lossfactor_mc", mc_samplesize, ".csv"))
        print("Lossfactor directory w/ required sample size already exists. Skipping...")
        return "Done"
    end

    # fix the seed and run the MC analysis
    Random.seed!(seed_used)
    results = runsim(model,
                     mc_draws_mutated,
                     tp_used == "all" || occursin("ISM", tp_used), # ism_used
                     omhversion != false && (tp_used == "all" || occursin("OMH", tp_used)), # omh_used
                     amocversion != false && (tp_used == "all" || occursin("AMOC", tp_used)), # amoc_used
                     safversion != false, # saf_used (always true)
                     ifelse(tp_used == "all" || occursin("AMAZ", tp_used), "Distribution", "none"), # AMAZ
                     ifelse(tp_used == "all" || occursin("WAIS", tp_used), "Distribution", "none"), # WAIS
                     save_lossfactor_as_csv = true, # save out the GDP loss as a CSV file
                     dir_lossfactor_csv = subfolder_path, # specify the location where CSV files are stored
                     country_iso3_subset = countries_with_marketcap # only save out GDP loss factors for countries in MSCI indices
                     ) # save all random variables)

    # save out MC results as a pickle-like HDF object
    if !isdir(string(dir_output, "/jld2_files"))
        mkdir(string(dir_output, "/jld2_files"))
    end

    save_object(string(dir_output, "/jld2_files/", parameter_id, ".jld2"),
                results)

    # return a Boolean indicating success
    return true
end


##############################################################
################## RUN SPECIFICATIONS OF INTEREST ############
##############################################################

# quick code to run the main specification (COMMENTED OUT)
# run_equityvaluation_mc(mc_samplesize_global)
# run_equityvaluation_mc(mc_samplesize_global, tp_used="none")

# Step 1: Create a DataFrame with all specifications to be run
# First we loop through the RCP-SSP pairs (first SSP2, then SSP5, all TPs) (12x)
# Then, we continue through the individual tipping points (SSP2-4.5) (6x)
#    NOTE: for OMH (which are deactivated in default spec), we need to switch omh from 'none' to 'default', otherwise this is equivalent to no TPs
# Then, different damage function specifications (SSP2-4.5, all TPs) (4x)
# Then, we run specs with and without tipping points a la Dietz et al 2021 (2x)
# Then, Cai et al. AMOC collapse (SSP2-4.5 and SSP2-8.5, all TPs) (2x)
specifications = DataFrame(
    rcp = vcat(repeat(["RCP4.5", "RCP3-PD/2.6", "RCP8.5"], inner=2, outer=2),
               repeat(["RCP4.5"], inner = 13),
               ["RCP8.5"]),
    ssp = vcat(repeat(["SSP2", "SSP5"], inner=6),
               repeat(["SSP2"], inner = 14)),
    tp_used = vcat(repeat(["all", "none"], outer=6), ["OMH", "PCF", "GIS", "WAIS", "AMAZ", "ISM"], repeat(["all", "none"], outer=3), ["all", "all"]),
    tdamage = vcat(repeat(["coacch_distribution"], inner=18), repeat(["bhm_distribution"], inner=4), repeat(["coacch_distribution"], inner=4)),
    persist = vcat(repeat([1.], inner=18), [0., 0., "Distribution", "Distribution"], repeat([1.], inner = 4)),
    omh = vcat(repeat(["none"], inner=12), ["default"], repeat(["none"], inner = 9), ["default"], repeat(["none"], inner = 3)),
    amoc = vcat(repeat(["none"], inner=22), ["IPSL", "none", "Cai", "Cai"]),
    saf = vcat(repeat(["Distribution mean"], inner=23), ["none", "Distribution mean", "Distribution mean"])
)

# export as a CSV file for inspection and reproducibility
CSV.write(string(dir_output_global, "/specifications_META.csv"), specifications)

# initiate an empty vector in which we store specifications that throw errors
specs_with_errors = Vector{Int}()

# lock to avoid data-racing
mylock = ReentrantLock()

# parallelize using @threads
# NOTE: as of August 2024, we do not use parallelization due to issues with @threads. UNCOMMENT the line below and COMMENT OUT the normal for-loop for parallelization
#Threads.@threads for ii in 1:nrow(specifications)
for ii in 1:nrow(specifications)
    try
        spec = specifications[ii, :]
        run_equityvaluation_mc(rcp = spec.rcp,
                               ssp = spec.ssp,
                               tp_used = spec.tp_used,
                               tdamage = spec.tdamage,
                               persist = spec.persist,
                               omh = spec.omh,
                               amoc = spec.amoc,
                               saf = spec.saf)
    catch e
        println("Error for specification ", ii, ": ", e)

        # add the specification to the vector using the reentrant lock
        lock(mylock) do
            push!(specs_with_errors, ii)
        end
    end
end

# report specifications with errors
println("Following specs were aborted due to errors: ", specs_with_errors, ". Please re-run single-threaded and/or inspect")

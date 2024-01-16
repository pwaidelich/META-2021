include("../lib/gdppc.jl")
include("../lib/damages.jl")
include("../lib/saverate.jl")

@defcomp Consumption begin
    # Variables
    gdppc_region = Variable(index=[time, region], unit="2010 USD PPP")
    gdppc_ratio_region = Variable(index=[time, region])
    gdppc_growth_region = Variable(index=[time, region])

    gdppc_growth = Variable(index=[time, country])
    gdppc = Variable(index=[time, country], unit="2010 USD PPP")
    conspc_preadj = Variable(index=[time, country], unit="2010 USD PPP") # previous year's
    conspc = Variable(index=[time, country], unit="2010 USD PPP")
    baseline_consumption_percap_percountry = Variable(index = [time, country], unit = "2010 USD PPP") # Counterfactual consumption per cap per country from SSPs

    # Parameters
    ssp = Parameter{String}()

    # Based on SSP
    convergerate_gdppc = Parameter()
    decayrate_gdppc = Parameter()
    popweights_region = Parameter(index=[region])

    # Based on historical data
    gdppc_2009 = Parameter(index=[country], unit="2010 USD PPP")
    saverate = Parameter(index=[country])

    # Based on damage specification
    seeds = Parameter{Int64}(index=[country])
    beta1 = Parameter(index=[country], unit="1/degC")
    beta2 = Parameter(index=[country], unit="1/degC^2")
    tempdamage = Variable(index=[time, country], unit = "%GDP")
    # parameters to implement BHM global function
    use_bhm_distribution = Parameter{Int64}()
    beta1_bhm = Parameter(unit="1/degC", default = 0.012718353)
    beta2_bhm = Parameter(unit="1/degC^2", default = -0.00048709)
    tempdamage_bhm = Variable(index=[time, country], unit = "%GDP")
    # parameters and variables to implement Waidelich et al 2024
    use_waid_distribution = Parameter{Int64}()
    beta0_waid = Parameter(index=[country])
    beta1_waid = Parameter(index=[country], unit="1/degC")
    beta2_waid = Parameter(index=[country], unit="1/degC^2")
    tempdamage_waid = Variable(index=[time, country], unit = "%GDP")
    T_AT = Parameter(index=[time], unit="degC")
    # parameters and variables to implement COACCH (van der Wijst et al 2023)
    beta1_coacch = Parameter(index=[country], unit="1/degC")
    beta2_coacch = Parameter(index=[country], unit="1/degC^2")
    tempdamage_coacch = Variable(index=[time, country], unit = "%GDP")
    # parameter determining which market damages are used (0 = Dietz, 1 = BHM, 2 = Waid et al, 3 = COACCH)
    tempdamage_switch = Parameter(default = 0)

    T_country_1990 = Parameter(index=[country], unit="degC")

    slrdamageconfig = Parameter{String}()
    slruniforms = Parameter(index=[country]) # only used for MC mode
    slrcoeff = Parameter(index=[country], unit="1/m")

    # Damage configuration
    damagepersist = Parameter(default=0.5)
    min_conspc = Parameter(unit="2010 USD PPP", default=1)
    extradamage = Parameter(index=[time, country])

    # Inputs from other components
    T_country = Parameter(index=[time, country], unit="degC")
    SLR = Parameter(index=[time], unit="m")

    function init(pp, vv, dd)
        isos = dim_keys(model, :country)

        for cc in dd.country
            if pp.seeds[cc] != 0
                betaboths = getbhmbetas(isos[cc], "distribution", pp.seeds[cc])
                pp.beta1[cc] = betaboths[1]
                pp.beta2[cc] = betaboths[2]
            end

            if pp.use_bhm_distribution != 0
                betaboths_bhm = getglobalbhmbetas_distribution()
                pp.beta1_bhm = betaboths_bhm[1]
                pp.beta2_bhm = betaboths_bhm[2]
            end

            if pp.use_waid_distribution != 0
                betaboths_waid = getwaidbetas(isos[cc], "waid_distribution")
                pp.beta0_waid[cc] = betaboths_waid[1]
                pp.beta1_waid[cc] = betaboths_waid[2]
                pp.beta2_waid[cc] = betaboths_waid[3]
            end
        end

        if pp.slrdamageconfig == "distribution"
            pp.slrcoeff = [getslrcoeff_distribution(isos[cc], slrdamage, pp.slruniforms[cc]) for cc in 1:length(isos)]
        end
    end

    function run_timestep(pp, vv, dd, tt)
        isos = dim_keys(model, :country)

        if is_first(tt)
            for rr in dd.region
                vv.gdppc_region[tt, rr] = getgdppc_ssp(dim_keys(model, :region)[rr], pp.ssp, 2010)
                if ismissing(vv.gdppc_region[tt, rr])
                    vv.gdppc_region[tt, rr] = 0
                end
                vv.gdppc_ratio_region[tt, rr] = 1
                vv.gdppc_growth_region[tt, rr] = 0 # ignored
            end

            for cc in dd.country
                gdppc = getgdppc(isos[cc], 2010)
                if ismissing(gdppc)
                    if isos[cc] == "DJI"
                        vv.gdppc[tt, cc] = getgdppc(isos[cc], 2011)
                        vv.gdppc_growth[tt, cc] = vv.gdppc[tt, cc] / pp.gdppc_2009[cc] - 1
                    else
                        vv.gdppc[tt, cc] = vv.gdppc_growth[tt, cc] = 0
                    end
                else
                    vv.gdppc[tt, cc] = gdppc
                    vv.gdppc_growth[tt, cc] = vv.gdppc[tt, cc] / pp.gdppc_2009[cc] - 1
                end

                vv.conspc_preadj[tt, cc] = (1-pp.saverate[cc])*pp.gdppc_2009[cc]
            end
        else
            for rr in dd.region
                vv.gdppc_region[tt, rr] = getgdppc_ssp(dim_keys(model, :region)[rr], pp.ssp, gettime(tt))
                vv.gdppc_ratio_region[tt, rr] = vv.gdppc_region[tt, rr] / vv.gdppc_region[TimestepIndex(1), rr]
                if gettime(tt) <= 2100
                    vv.gdppc_growth_region[tt, rr] = vv.gdppc_ratio_region[tt, rr] / vv.gdppc_ratio_region[tt-1, rr] - 1
                else
                    vv.gdppc_growth_region[tt, rr] = (1-pp.convergerate_gdppc-pp.decayrate_gdppc)*vv.gdppc_growth_region[tt-1, rr]+pp.decayrate_gdppc*sum(vv.gdppc_growth_region[tt-1, :] .* pp.popweights_region)
                end
            end

            for cc in dd.country
                region = getregion(isos[cc])
                if ismissing(region)
                    growth = 0
                else
                    rr = findfirst(dim_keys(model, :region) .== region)
                    growth = vv.gdppc_growth_region[tt, rr]
                end

                vv.gdppc_growth[tt, cc] = growth
                vv.gdppc[tt, cc] = vv.gdppc[tt-1, cc] * (1 + growth)

                if vv.conspc[tt-1, cc] == 0 || vv.gdppc[tt, cc] == 0
                    vv.conspc_preadj[tt, cc] = 0
                else
                    vv.conspc_preadj[tt, cc] = pp.damagepersist*(1-pp.saverate[cc])*vv.gdppc[tt-1, cc]+(1-pp.damagepersist)*vv.conspc[tt-1, cc]
                end
            end
        end

        for cc in dd.country
            # NOTE: signs are reversed because for Dietz et al 2021 country-specific damage function, positive value = damage
            vv.tempdamage[tt, cc] = - pp.beta1[cc]*(pp.T_country[tt, cc]-pp.T_country_1990[cc]) - pp.beta2[cc]*(pp.T_country[tt, cc]-pp.T_country_1990[cc])^2
            vv.tempdamage_bhm[tt, cc] = pp.beta1_bhm*(pp.T_country[tt, cc]-pp.T_country_1990[cc]) + pp.beta2_bhm*(pp.T_country[tt, cc]^2-pp.T_country_1990[cc]^2)
            vv.tempdamage_waid[tt, cc] = pp.beta0_waid[cc] + pp.beta1_waid[cc]*pp.T_AT[tt] + pp.beta2_waid[cc]*(pp.T_AT[tt]^2)
            # we use 2010 global warming as baseline - van der Wijst et al 2023 use 1986-2005 average (see their Fig. 1) XX
            vv.tempdamage_coacch[tt, cc] = (1/100) * (- pp.beta1_coacch[cc]*(pp.T_AT[tt] - pp.T_AT[TimestepIndex(1)]) - pp.beta2_coacch[cc]*((pp.T_AT[tt] - pp.T_AT[TimestepIndex(1)])^2))

            if pp.tempdamage_switch == 0
                vv.conspc[tt, cc] = vv.conspc_preadj[tt, cc]*(1+vv.gdppc_growth[tt, cc]+vv.tempdamage[tt, cc])*(1-pp.SLR[tt]*pp.slrcoeff[cc])*(1 - pp.extradamage[tt, cc])
            elseif pp.tempdamage_switch == 1
                vv.conspc[tt, cc] = vv.conspc_preadj[tt, cc]*(1+vv.gdppc_growth[tt, cc]+vv.tempdamage_bhm[tt, cc])*(1-pp.SLR[tt]*pp.slrcoeff[cc])*(1 - pp.extradamage[tt, cc])
            elseif pp.tempdamage_switch == 2
                vv.conspc[tt, cc] = vv.conspc_preadj[tt, cc]*(1+vv.gdppc_growth[tt, cc])*(1+vv.tempdamage_waid[tt, cc])*(1-pp.SLR[tt]*pp.slrcoeff[cc])*(1 - pp.extradamage[tt, cc])
            elseif pp.tempdamage_switch == 3
                vv.conspc[tt, cc] = vv.conspc_preadj[tt, cc]*(1+vv.gdppc_growth[tt, cc])*(1+vv.tempdamage_coacch[tt, cc])*(1-pp.SLR[tt]*pp.slrcoeff[cc])*(1 - pp.extradamage[tt, cc])
            end

            # Compute baseline consumption per capita without damages
            vv.baseline_consumption_percap_percountry[tt,cc] = (1-pp.saverate[cc])*vv.gdppc[tt, cc]

            if vv.conspc[tt, cc] <= pp.min_conspc && vv.conspc[TimestepIndex(1), cc] != 0
                vv.conspc[tt, cc] = pp.min_conspc
            end
        end
    end
end

function addConsumption(model, tdamage, slrdamage, ssp)
    if tdamage ∉ ["none", "distribution", "pointestimate", "low", "high",
                  "waid_pointestimate", "waid_distribution",
                  "bhm_pointestimate", "bhm_distribution",
                  "coacch_central"]
        throw(ArgumentError("Unknown Consumption tdamage"))
    end
    if slrdamage ∉ ["none", "distribution", "mode", "low", "high"]
        throw(ArgumentError("Unknown Consumption slrdamage"))
    end

    cons = add_comp!(model, Consumption)

    cons[:ssp] = ssp

    sspextend = CSV.read("../data/sspextend-gdppc.csv", DataFrame)
    cons[:convergerate_gdppc] = sspextend[sspextend.SSP .== ssp, "Convergence rate"][1]
    cons[:decayrate_gdppc] = sspextend[sspextend.SSP .== ssp, "Decay rate"][1]

    cons[:popweights_region] = [getpopweight_ssp(region, ssp) for region in dim_keys(model, :region)]

    isos = dim_keys(model, :country)

    # unless tdamage is "none" or "distribution", we use Dietz et al point estimates
    if tdamage == "none"
        cons[:seeds] = [0 for iso in isos]
        cons[:beta1] = zeros(length(isos))
        cons[:beta2] = zeros(length(isos))
    elseif tdamage != "distribution"
        betaboths = [getbhmbetas(iso, tdamage) for iso in isos]
        cons[:seeds] = [0 for iso in isos]
        cons[:beta1] = [betaboth[1] for betaboth in betaboths]
        cons[:beta2] = [betaboth[2] for betaboth in betaboths]
    else
        cons[:beta1] = zeros(length(isos))
        cons[:beta2] = zeros(length(isos))
    end

    # if BHM distribution is set, we draw - otherwise, we use point estimates (= default parameter values)
    if tdamage == "bhm_distribution"
        cons[:beta1_bhm] = 0.
        cons[:beta2_bhm] = 0.
        cons[:use_bhm_distribution] = 1
    else
        cons[:use_bhm_distribution] = 0
    end

    # if Waidelich et al distribution is set, we draw - otherwise, we populate with point estimates
    if tdamage == "waid_distribution"
        cons[:use_waid_distribution] = 1
        cons[:beta0_waid] = zeros(length(isos))
        cons[:beta1_waid] = zeros(length(isos))
        cons[:beta2_waid] = zeros(length(isos))
    else
        cons[:use_waid_distribution] = 0
        betaboths_waid = [getwaidbetas(iso, "waid_pointestimate") for iso in isos]
        cons[:beta0_waid] = [betaboth[1] for betaboth in betaboths_waid]
        cons[:beta1_waid] = [betaboth[2] for betaboth in betaboths_waid]
        cons[:beta2_waid] = [betaboth[3] for betaboth in betaboths_waid]
    end

    # COACCH betas default to central value (= 50th quantile)
    betaboths_coacch = [getcoacchbetas(iso) for iso in isos]
    cons[:beta1_coacch] = [betaboth[1] for betaboth in betaboths_coacch]
    cons[:beta2_coacch] = [betaboth[2] for betaboth in betaboths_coacch]

    if slrdamage == "none"
        cons[:slrcoeff] = zeros(length(isos))
        cons[:slruniforms] = zeros(length(isos))
    elseif slrdamage != "distribution"
        cons[:slrcoeff] = [getslrcoeff(iso, slrdamage) for iso in isos]
        cons[:slruniforms] = zeros(length(isos))
    else
        cons[:slrcoeff] = zeros(length(isos))  # to be filled by init
        cons[:slruniforms] = zeros(length(isos))  # to be filled by monte carlo
    end

    cons[:slrdamageconfig] = slrdamage

    cons[:saverate] = [getsaverate(iso) for iso in isos]
    gdppc_2009 = [getgdppc(iso, 2009) for iso in isos]
    gdppc_2009[ismissing.(gdppc_2009)] .= 0
    gdppc_2009[isos .== "DJI"] .= 2700
    cons[:gdppc_2009] = gdppc_2009
    cons[:T_country_1990] = [gettemp1990(iso) for iso in isos]
    cons[:extradamage] = zeros(dim_count(model, :time), dim_count(model, :country))

    if tdamage ∈ ["bhm_pointestimate", "bhm_distribution"]
        cons[:tempdamage_switch] = 1
    elseif tdamage ∈ ["waid_pointestimate", "waid_distribution"]
        cons[:tempdamage_switch] = 2
    elseif tdamage ∈ ["coacch_central"]
        cons[:tempdamage_switch] = 3
    end

    cons
end

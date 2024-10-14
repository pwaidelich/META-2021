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
    lossfactor_conspc = Variable(index=[time, country])
    lossfactor_temp = Variable(index=[time, country])
    lossfactor_persistence = Variable(index=[time, country])
    lossfactor_SLR = Variable(index=[time, country])
    lossfactor_extradamages = Variable(index=[time, country])

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
    #use_bhm_distribution = Parameter{Int64}()
    beta1_bhm = Parameter(unit="1/degC", default = 0.012718353)
    beta2_bhm = Parameter(unit="1/degC^2", default = -0.00048709)
    tempdamage_bhm = Variable(index=[time, country], unit = "%GDP")
    # parameters and variables to implement Waidelich et al 2024
    use_waid_distribution = Parameter{Int64}()
    #beta0_waid = Parameter(index=[country])
    beta1_waid = Parameter(index=[country], unit="1/degC")
    beta2_waid = Parameter(index=[country], unit="1/degC^2")
    tempdamage_waid = Variable(index=[time, country], unit = "%GDP")
    T_AT = Parameter(index=[time], unit="degC")
    T_AT_baseline_waid = Parameter(unit="degC", default = 0.84)
    #T_AT_cap_waid = Parameter(unit="degC", default = 6.181668)
    # parameters and variables to implement COACCH (van der Wijst et al 2023)
    use_coacch_distribution = Parameter{Int64}()
    seed_coacch = Parameter{Int64}(default = 0)
    beta1_coacch = Parameter(index=[country], unit="1/degC")
    beta2_coacch = Parameter(index=[country], unit="1/degC^2")
    alpha_coacch = Parameter(index=[country])
    T_AT_1986_2005 = Parameter(unit="degC", default = 0.776462745) # calculated via Berkeley Earth Annual Anomalies (Land + Ocean anomaly using air temperature above sea ice)
    tempdamage_coacch = Variable(index=[time, country], unit = "%GDP")
    deltaTglobal_country_AMOC = Parameter(index=[time, country], unit="degC")
    # Kahn damages
    beta_plus_kahn = Parameter(unit="1/degC", default = -0.058580645) # positive deviation from historical norm
    beta_minus_kahn = Parameter(unit="1/degC", default = -0.052) # negative deviation from historical norm
    T_country_1981_2010 = Variable(index =[country], unit = "degC")
    T_country_kahn_norm = Variable(index=[time, country], unit = "degC")
    T_country_kahn_deviation = Variable(index = [time, country], unit = "degC")
    kahn_norm_end = Variable(index = [time])
    tempdamage_kahn = Variable(index=[time, country], unit = "%GDP")

    # parameter determining which market damages are used (0 = Dietz, 1 = BHM, 2 = Waid et al, 3 = COACCH)
    tempdamage_switch = Parameter(default = 1)

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

    # AMOC-based SLR adder XX
    #SLR_AMOC_adder = Parameter(index = [time], unit = "m")
    extradamage_AMOC = Parameter(index=[time, country])

    function init(pp, vv, dd)
        isos = dim_keys(model, :country)

        #if pp.use_bhm_distribution != 0
        #    betaboths_bhm = getglobalbhmbetas_distribution()
        #    pp.beta1_bhm = betaboths_bhm[1]
        #    pp.beta2_bhm = betaboths_bhm[2]
        #end

        for cc in dd.country
            if pp.seeds[cc] != 0
                betaboths = getbhmbetas(isos[cc], "distribution", pp.seeds[cc])
                pp.beta1[cc] = betaboths[1]
                pp.beta2[cc] = betaboths[2]
            end

            if pp.use_waid_distribution != 0
                betaboths_waid = getwaidbetas(isos[cc], "waid_distribution")
                #pp.beta0_waid[cc] = betaboths_waid[1]
                pp.beta1_waid[cc] = betaboths_waid[1]
                pp.beta2_waid[cc] = betaboths_waid[2]
            end

            if pp.use_coacch_distribution == 1
                pp.alpha_coacch[cc] = getcoacchalpha_distribution(isos[cc], pp.seed_coacch)
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

                # calculate the 2001-2009 temperatures by linear interpolation
                vv.T_country_kahn_norm[tt, cc] = (2/3)*pp.T_country_1990[cc] + 1/3*pp.T_country[tt, cc]
                vv.T_country_kahn_deviation[tt, cc] = pp.T_country[tt, cc] - vv.T_country_kahn_norm[tt, cc]

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

                # if we are prior to 2039 (timestep 30), we calculate 30-year average based on 1981-2010 average
                if gettime(tt) <= 2040
                    # we add the previous year and make the average of the previous 29 years one increment higher (assuming linear increase)
                    #vv.kahn_norm_end[tt] = TimestepIndex(tt-1).t
                    # NOTE: below we access the integer of tt by getting the year and subtracting 2009 (1 for 2010, 2 for 2011 and so on). SOMEWHAT HACKISH XX
                    #vv.T_country_kahn_norm[tt, cc] = vv.T_country_kahn_norm[tt-1, cc]*(gettime(tt) - 2009 - 1)/(gettime(tt) - 2009) + pp.T_country[tt-1, cc]*(1/(gettime(tt) - 2009))
                #else
                    vv.T_country_kahn_norm[tt, cc] = vv.T_country_kahn_norm[tt-1, cc] - vv.T_country_kahn_norm[tt-1, cc]/30 + pp.T_country[tt-1, cc]/30
                #elseif gettime(tt) <= 2040
                #    vv.T_country_kahn_norm[tt, cc] = vv.T_country_kahn_norm[tt-1, cc] - pp.T_country[TimestepIndex(1), cc]/30 + pp.T_country[tt-1, cc]/30
                else
                    vv.T_country_kahn_norm[tt, cc] = vv.T_country_kahn_norm[tt-1, cc] - pp.T_country[tt-30, cc]/30 + pp.T_country[tt-1, cc]/30
                end
                vv.T_country_kahn_deviation[tt, cc] = pp.T_country[tt, cc] - vv.T_country_kahn_norm[tt, cc]

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
            vv.tempdamage_waid[tt, cc] = pp.beta1_waid[cc]*(pp.T_AT[tt] - pp.T_AT_baseline_waid) + pp.beta2_waid[cc]*(pp.T_AT[tt] - pp.T_AT_baseline_waid)^2
            # van der Wijst et al 2023 damage functions use global warming over 1986-2005 average and return results in %GDP (> 0 = damage)
            # NOTE: we floor warming at the baseline level because COACCH functions are not calibrated for negative warming
            vv.tempdamage_coacch[tt, cc] = (1/100) * pp.alpha_coacch[cc] * (- pp.beta1_coacch[cc]*max(0, pp.T_AT[tt] + pp.deltaTglobal_country_AMOC[tt, cc] - pp.T_AT_1986_2005) - pp.beta2_coacch[cc]*(max(0, pp.T_AT[tt] + pp.deltaTglobal_country_AMOC[tt, cc] - pp.T_AT_1986_2005)^2))
            vv.tempdamage_kahn[tt, cc] = vv.T_country_kahn_deviation[tt, cc] > 0 ? pp.beta_plus_kahn*vv.T_country_kahn_deviation[tt, cc] : pp.beta_minus_kahn*vv.T_country_kahn_deviation[tt, cc]

            if pp.tempdamage_switch == 0
                vv.conspc[tt, cc] = vv.conspc_preadj[tt, cc]*(1+vv.gdppc_growth[tt, cc]+vv.tempdamage[tt, cc])*(1-pp.SLR[tt]*pp.slrcoeff[cc])*(1 - pp.extradamage[tt, cc])*(1 - pp.extradamage_AMOC[tt, cc])
            elseif pp.tempdamage_switch == 1
                vv.conspc[tt, cc] = vv.conspc_preadj[tt, cc]*(1+vv.gdppc_growth[tt, cc]+vv.tempdamage_bhm[tt, cc])*(1-pp.SLR[tt]*pp.slrcoeff[cc])*(1 - pp.extradamage[tt, cc])*(1 - pp.extradamage_AMOC[tt, cc])
            elseif pp.tempdamage_switch == 2
                vv.conspc[tt, cc] = vv.conspc_preadj[tt, cc]*(1+vv.gdppc_growth[tt, cc])*(1+vv.tempdamage_waid[tt, cc])*(1-pp.SLR[tt]*pp.slrcoeff[cc])*(1 - pp.extradamage[tt, cc])*(1 - pp.extradamage_AMOC[tt, cc])
            elseif pp.tempdamage_switch == 3
                vv.conspc[tt, cc] = vv.conspc_preadj[tt, cc]*(1+vv.gdppc_growth[tt, cc])*(1+vv.tempdamage_coacch[tt, cc])*(1-pp.SLR[tt]*pp.slrcoeff[cc])*(1 - pp.extradamage[tt, cc])*(1 - pp.extradamage_AMOC[tt, cc])
            elseif pp.tempdamage_switch == 4
                vv.conspc[tt, cc] = vv.conspc_preadj[tt, cc]*(1+vv.gdppc_growth[tt, cc])*(1+vv.tempdamage_kahn[tt, cc])*(1-pp.SLR[tt]*pp.slrcoeff[cc])*(1 - pp.extradamage[tt, cc])*(1 - pp.extradamage_AMOC[tt, cc])
            end

            # Compute baseline consumption per capita without damages
            vv.baseline_consumption_percap_percountry[tt,cc] = (1-pp.saverate[cc])*vv.gdppc[tt, cc]

            if vv.conspc[tt, cc] <= pp.min_conspc && vv.conspc[TimestepIndex(1), cc] != 0
                vv.conspc[tt, cc] = pp.min_conspc
            end

            # Compute relative loss in consumption per capita compared to baseline
            vv.lossfactor_conspc[tt, cc] = vv.conspc[tt, cc]/vv.baseline_consumption_percap_percountry[tt,cc]
            vv.lossfactor_temp[tt, cc] = 1 - pp.beta1[cc]*(pp.T_country[tt, cc] - pp.T_country_1990[cc]) - pp.beta2[cc]*(pp.T_country[tt, cc] - pp.T_country_1990[cc])^2
            vv.lossfactor_persistence[tt,cc] = (vv.conspc_preadj[tt, cc]*(1+(vv.gdppc_growth[tt, cc])))/vv.baseline_consumption_percap_percountry[tt,cc]
            vv.lossfactor_SLR[tt,cc] = (1-pp.SLR[tt]*pp.slrcoeff[cc])
            vv.lossfactor_extradamages[tt,cc] = (1 - pp.extradamage[tt, cc])
        end
    end
end

function addConsumption(model, tdamage, slrdamage, ssp)
    if tdamage ∉ ["none", "distribution", "pointestimate", "low", "high",
                  "waid_pointestimate", "waid_distribution",
                  "bhm_pointestimate", "bhm_distribution",
                  "coacch_central", "coacch_distribution",
                  "kahn_pointestimate"]
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
    end

    # if Waidelich et al distribution is set, we draw - otherwise, we populate with point estimates
    if tdamage == "waid_distribution"
        cons[:use_waid_distribution] = 1
        #cons[:beta0_waid] = zeros(length(isos))
        cons[:beta1_waid] = zeros(length(isos))
        cons[:beta2_waid] = zeros(length(isos))
    else
        cons[:use_waid_distribution] = 0
        betaboths_waid = [getwaidbetas(iso, "waid_pointestimate") for iso in isos]
        #cons[:beta0_waid] = [betaboth[1] for betaboth in betaboths_waid]
        cons[:beta1_waid] = [betaboth[1] for betaboth in betaboths_waid]
        cons[:beta2_waid] = [betaboth[2] for betaboth in betaboths_waid]
    end

    # COACCH betas default to central value (= 50th quantile) unless
    betaboths_coacch = [getcoacchbetas(iso) for iso in isos]
    cons[:beta1_coacch] = [betaboth[1] for betaboth in betaboths_coacch]
    cons[:beta2_coacch] = [betaboth[2] for betaboth in betaboths_coacch]

    if tdamage == "coacch_distribution"
        cons[:use_coacch_distribution] = 1
        cons[:alpha_coacch] = zeros(length(isos))
    else
        cons[:use_coacch_distribution] = 0
        cons[:alpha_coacch] = ones(length(isos))
    end

    # if COACCH distribution is set, we overwrite with zeros and draw betas in the component


    cons[:deltaTglobal_country_AMOC] = zeros(dim_count(model, :time), dim_count(model, :country))

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

    # set the AMOC SLR adder to zero XX
    #cons[:SLR_AMOC_adder] = zeros(dim_count(model, :time))

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
    elseif tdamage ∈ ["coacch_central", "coacch_distribution"]
        cons[:tempdamage_switch] = 3
    elseif tdamage ∈ ["kahn_pointestimate"]
        cons[:tempdamage_switch] = 4
    end

    cons
end

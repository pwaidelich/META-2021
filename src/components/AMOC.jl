@defcomp AMOC begin
    country = Index()

    # Variables
    p_AMOC = Variable(index=[time], unit="fraction")
    I_AMOC = Variable{Bool}(index=[time])
    deltaT_country_AMOC = Variable(index=[time, country], unit="degC")
    T_country_AMOC = Variable(index=[time, country], unit="degC")

    # Parameters
    T_AT = Parameter(index=[time], unit="degC")
    f_AMOC = Parameter(index=[time])

    uniforms = Parameter(index=[time])
    b_AMOC = Parameter()
    Delta_AMOC = Parameter(unit="year", default=35)
    max_deltaT_country_AMOC = Parameter(index=[country], unit="degC")

    T_country_base = Parameter(index=[time, country], unit="degC")

    # new parameter to emulate AMOC impacts when using COACCH damage functions
    scale_country = Parameter(index=[time, country])
    deltaTglobal_country_AMOC = Variable(index=[time, country], unit="degC")
    use_AMOC_Cai_impacts = Parameter(default = 0.)
    gdp_impact_proportional = Parameter(index=[country])
    extradamage_AMOC = Variable(index=[time, country])
    Delta_AMOC_Cai = Parameter(unit="year", default=50)
    max_GDPimpact_Cai = Parameter(default=0.15)

    #SLR_AMOC_adder = Variable(index=[time], unit="m")

    function run_timestep(pp, vv, dd, tt)
        # if we use catastropic AMOC impacts, set hazard rate much lower (to ensure <=10% probability of tipping as per IPCC)
        if pp.use_AMOC_Cai_impacts == 1.
            # NOTE: we take hazard rate from Cai et al. (2016), which implies roughly a 5% chance of tipping within 21st century under RCP4.5 MC mean
            pp.b_AMOC = 0.00063064
            # if we use the Cai et al. AMOC collapse, then one degC is substracted from global warming over pre-industrial
            vv.p_AMOC[tt] = min((1-exp(-pp.b_AMOC*max(0, pp.T_AT[tt] - 1)))*pp.f_AMOC[tt], 1)
        elseif pp.use_AMOC_Cai_impacts == 0.
            vv.p_AMOC[tt] = min((1-exp(-pp.b_AMOC*pp.T_AT[tt]))*pp.f_AMOC[tt], 1)
        end



        if is_first(tt)
            vv.I_AMOC[tt] = pp.uniforms[tt] < vv.p_AMOC[tt]

            for cc in dd.country
                vv.deltaT_country_AMOC[tt, cc] = vv.I_AMOC[tt] ? pp.max_deltaT_country_AMOC[cc] / pp.Delta_AMOC : 0

                # NOTE: for Cai et al. 2016 impacts of AMOC collapse, we use different transition period and scale impacts based on maximum %GDP loss
                vv.extradamage_AMOC[tt, cc] = vv.I_AMOC[tt] ? (pp.gdp_impact_proportional[cc] * pp.max_GDPimpact_Cai) / pp.Delta_AMOC_Cai : 0
            end
        else
            vv.I_AMOC[tt] = vv.I_AMOC[tt-1] || (pp.uniforms[tt] < vv.p_AMOC[tt])

            for cc in dd.country
                if vv.I_AMOC[tt] && abs(vv.deltaT_country_AMOC[tt-1, cc]) < abs(pp.max_deltaT_country_AMOC[cc])
                    vv.deltaT_country_AMOC[tt, cc] = vv.deltaT_country_AMOC[tt-1, cc] + pp.max_deltaT_country_AMOC[cc] / pp.Delta_AMOC
                    if abs(vv.deltaT_country_AMOC[tt, cc]) > abs(pp.max_deltaT_country_AMOC[cc])
                        vv.deltaT_country_AMOC[tt, cc] = pp.max_deltaT_country_AMOC[cc]
                    end

                    # Cai et al. 2016 AMOC collapse impacts
                    vv.extradamage_AMOC[tt, cc] = vv.extradamage_AMOC[tt-1, cc] + (pp.gdp_impact_proportional[cc] * pp.max_GDPimpact_Cai) / pp.Delta_AMOC_Cai
                    if abs(vv.extradamage_AMOC[tt, cc]) > abs(pp.gdp_impact_proportional[cc] * pp.max_GDPimpact_Cai)
                        vv.extradamage_AMOC[tt, cc] = pp.gdp_impact_proportional[cc] * pp.max_GDPimpact_Cai
                    end

                else
                    vv.deltaT_country_AMOC[tt, cc] = vv.deltaT_country_AMOC[tt-1, cc]

                    vv.extradamage_AMOC[tt, cc] = vv.extradamage_AMOC[tt-1, cc]
                end
            end
        end

        for cc in dd.country
            vv.T_country_AMOC[tt, cc] = pp.T_country_base[tt, cc] + vv.deltaT_country_AMOC[tt, cc]

            # emulate what global temperature would cause the same shift in country cc's national temperature
            if pp.use_AMOC_Cai_impacts == 0.
                vv.deltaTglobal_country_AMOC[tt, cc] = vv.deltaT_country_AMOC[tt, cc]/pp.scale_country[tt, cc]
                vv.extradamage_AMOC[tt, cc] = 0.
            else
                # if we use Cai et al 2016 AMOC collapse impacts, simply set the DeltaTglobal to zero
                vv.deltaTglobal_country_AMOC[tt, cc] = 0.
            end
        end
    end
end

b_AMOC_calibs = Dict{String, Float64}("Hadley" => 0.135, "BCM" => 0.611,
	                              "IPSL" => 0.54, "HADCM" => 1.6)

function addAMOC(model, calibration; before=nothing, after=nothing)
    if calibration ∉ keys(b_AMOC_calibs)
        throw(ArgumentError("Unknown AMOC model calibration"))
    end

    params = CSV.read("../data/AMOCparams.csv", DataFrame)
    params_GDPimpacts = CSV.read("../data/AMOCparams_GDPimpacts.csv", DataFrame)

    amoc = add_comp!(model, AMOC, before=before, after=after)
    amoc[:b_AMOC] = b_AMOC_calibs[calibration]
    amoc[:max_deltaT_country_AMOC] = params[!, calibration]

    amoc[:gdp_impact_proportional] = params_GDPimpacts[!, "gdp_impact_proportional"]

    amoc[:f_AMOC] = ones(dim_count(model, :time))

    amoc
end

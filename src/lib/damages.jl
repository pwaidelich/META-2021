import Random

bhmbetas = CSV.read("../data/BHMbetas.csv", DataFrame)
amocparams = CSV.read("../data/AMOCparams.csv", DataFrame)
waidbetas = CSV.read("../data/WaidelichEtAlbetas.csv", DataFrame)
coacchbetas = CSV.read("../data/COACCHbetas.csv", DataFrame)

function getbhmbetas(iso, option, seed=nothing)
    if option == "pointestimate"
        beta1 = bhmbetas.beta1[bhmbetas.iso .== iso][1]
        beta2 = bhmbetas.beta2[bhmbetas.iso .== iso][1]
    elseif option == "low"
        beta1 = bhmbetas.beta1[bhmbetas.iso .== iso][1] - 1.96*sqrt(bhmbetas.var11[bhmbetas.iso .== iso][1])
        beta2 = bhmbetas.beta2[bhmbetas.iso .== iso][1] - 1.96*sqrt(bhmbetas.var22[bhmbetas.iso .== iso][1])
    elseif option == "high"
        beta1 = bhmbetas.beta1[bhmbetas.iso .== iso][1] + 1.96*sqrt(bhmbetas.var11[bhmbetas.iso .== iso][1])
        beta2 = bhmbetas.beta2[bhmbetas.iso .== iso][1] + 1.96*sqrt(bhmbetas.var22[bhmbetas.iso .== iso][1])
    else
        if seed != nothing
            Random.seed!(seed)
        end
        mu1, mu2 = getbhmbetas(iso, "pointestimate")
        var11 = bhmbetas.var11[bhmbetas.iso .== iso][1]
        var12 = bhmbetas.var12[bhmbetas.iso .== iso][1]
        var22 = bhmbetas.var22[bhmbetas.iso .== iso][1]
        mvn = MvNormal([mu1, mu2], [var11 var12; var12 var22])
        beta1, beta2 = rand(mvn)
    end

    return beta1, beta2
end

function getglobalbhmbetas_distribution(seed=nothing)
    if seed != nothing
        Random.seed!(seed)
    end
    mu1, mu2 = (0.012718353, -0.00048709)
    var11 = 0.0000143458
    var12 = -0.000000375824
    var22 = 0.0000000140157
    mvn = MvNormal([mu1, mu2], [var11 var12; var12 var22])
    beta1_bhm, beta2_bhm = rand(mvn)
    return beta1_bhm, beta2_bhm
end

function getwaidbetas(iso, option, seed=nothing)
    if option == "waid_pointestimate"
        beta1_waid = waidbetas.beta1_waid[waidbetas.iso .== iso][1]
        beta2_waid = waidbetas.beta2_waid[waidbetas.iso .== iso][1]
    elseif option == "waid_distribution"
        if seed != nothing
            Random.seed!(seed)
        end

        mu1, mu2 = getwaidbetas(iso, "waid_pointestimate")
        var11 = waidbetas.var11[waidbetas.iso .== iso][1]

        if mu2 != 0.
            var12 = waidbetas.var12[waidbetas.iso .== iso][1]
            var22 = waidbetas.var22[waidbetas.iso .== iso][1]
            mvn = MvNormal([mu1, mu2], [var11 var12; var12 var22])
            beta1_waid, beta2_waid = rand(mvn)
        else
            uvn = Normal(mu1, var11)
            beta1_waid = rand(uvn)
            beta2_waid = 0.
        end
    end

    return beta1_waid, beta2_waid
end

function getcoacchbetas(iso)
    beta1_coacch = coacchbetas.no_slr_b1[coacchbetas.iso .== iso][1]
    beta2_coacch = coacchbetas.no_slr_b2[coacchbetas.iso .== iso][1]

    return beta1_coacch, beta2_coacch
end

function gettemp1990(iso)
    return amocparams[amocparams."Country code" .== iso, "1990 temp"][1]
end

slrcoeffs = CSV.read("../data/SLRcoeffs.csv", DataFrame)

function getslrcoeff(iso, option::String)
    if option == "mode"
        return slrcoeffs.mode[slrcoeffs.ISO .== iso][1]
    elseif option == "low"
        return slrcoeffs.low[slrcoeffs.ISO .== iso][1]
    elseif option == "high"
        return slrcoeffs.hig[slrcoeffs.ISO .== iso][1]
    else
        throw(ArgumentError("Distribution implemented in getslrcoeff_distribution"))
    end
end

function getslrcoeff_distribution(iso, option::String, qq)
    low = slrcoeffs.low[slrcoeffs.ISO .== iso][1]
    mode = slrcoeffs.mode[slrcoeffs.ISO .== iso][1]
    high = slrcoeffs.hig[slrcoeffs.ISO .== iso][1]
    quantile(TriangularDist(low, high, mode), qq)
end

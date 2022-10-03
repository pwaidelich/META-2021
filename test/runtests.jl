using Test

@testset "MimiMETA-all" begin
    include("test_RCP.jl")
    include("test_CO2Model.jl")
    include("test_CH4Model.jl")
    include("test_SAF.jl")
    include("test_TemperatureModel.jl")
    include("test_PCF.jl")
    include("test_OMH.jl")
    include("test_AmazonDieback.jl")
    include("test_GISModel.jl")
    include("test_WAIS.jl")
    include("test_SLR.jl")
    include("test_ISM.jl")
    include("test_AMOC.jl")
    include("test_PatternScaling.jl")
    include("test_Interactions.jl")
    include("test_Consumption.jl")
    include("test_NonMarketDamages.jl")
    include("test_Utility.jl")
end

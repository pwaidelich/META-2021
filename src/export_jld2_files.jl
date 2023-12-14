using CSV
using DataFrames
using JLD2

# Replace "your_file.jld2" with the path to your .jld2 file
# Replace "your_dataframe" with the name of the variable storing the DataFrame in the .jld2 file

using JLD2

years = collect(2010:2200)

# map the directory
dir_meta_results = "C:/Users/pwaidelich/Downloads/GitHub - Local/equity_var_climate/data/META MC results"

# Open the .jld2 file
file = jldopen(string(dir_meta_results, "/jld2_files/500_RCP45SSP2_persistDistribution_tipall_tdamageBHM.jld2"), "r")  # Replace "your_file.jld2" with your file's name

# inspect the keys in it
keys(file)
# -> "single_stored_object"


# NOTE: the resulting object is a vector with length = # of MC draws. Each element then is a dictionary

function extract_global_temperature(filepath = String[])

    # load the results
    @load filepath single_stored_object

    df = DataFrame(mc_run = Int64[], year = Int64[], T_AT = Float64[])

    for jj in 1:length(single_stored_object)

        df_to_add = DataFrame(mc_run = jj, year = years, T_AT = single_stored_object[jj][:TemperatureModel_T_AT])

        df = vcat(df, df_to_add)
    end

    return df
end

# RCP4.5 with TPs
df_T_AT_tipall_RCP45 = extract_global_temperature(string(dir_meta_results, "/jld2_files/500_RCP45SSP2_persistDistribution_tipall_tdamageBHM.jld2"))
CSV.write(string(dir_meta_results, "/T_AT/T_AT_n500_RCP45SSP2_tipall.csv"), df_T_AT_tipall_RCP45)

# RCP 4.5 w/o TPs
df_T_AT_tipnone_RCP45 = extract_global_temperature(string(dir_meta_results, "/jld2_files/500_RCP45SSP2_persistDistribution_tipnone_tdamageBHM.jld2"))
CSV.write(string(dir_meta_results, "/T_AT/T_AT_n500_RCP45SSP2_tipnone.csv"), df_T_AT_tipnone_RCP45)


# RCP8.5 with TPs
df_T_AT_tipall_RCP85 = extract_global_temperature(string(dir_meta_results, "/jld2_files/500_RCP85SSP2_persistDistribution_tipall_tdamageBHM.jld2"))
CSV.write(string(dir_meta_results, "/T_AT/T_AT_n500_RCP85SSP2_tipall.csv"), df_T_AT_tipall_RCP85)

# RCP 8.5 w/o TPs
df_T_AT_tipnone_RCP85 = extract_global_temperature(string(dir_meta_results, "/jld2_files/500_RCP85SSP2_persistDistribution_tipnone_tdamageBHM.jld2"))
CSV.write(string(dir_meta_results, "/T_AT/T_AT_n500_RCP85SSP2_tipnone.csv"), df_T_AT_tipnone_RCP85)

# RCP2.6 with TPs
df_T_AT_tipall_RCP26 = extract_global_temperature(string(dir_meta_results, "/jld2_files/500_RCP3-PD26SSP2_persistDistribution_tipall_tdamageBHM.jld2"))
CSV.write(string(dir_meta_results, "/T_AT/T_AT_n500_RCP3-PD26SSP2_tipall.csv"), df_T_AT_tipall_RCP26)

# RCP 2.6 w/o TPs
df_T_AT_tipnone_RCP26 = extract_global_temperature(string(dir_meta_results, "/jld2_files/500_RCP3-PD26SSP2_persistDistribution_tipnone_tdamageBHM.jld2"))
CSV.write(string(dir_meta_results, "/T_AT/T_AT_n500_RCP3-PD26SSP2_tipnone.csv"), df_T_AT_tipnone_RCP26)

# repeat the same for the global sea-level rise variable
function extract_global_SLR(filepath = String[])

    # load the results
    @load filepath single_stored_object

    df = DataFrame(mc_run = Int64[], year = Int64[], SLR = Float64[])

    for jj in 1:length(single_stored_object)

        df_to_add = DataFrame(mc_run = jj, year = years, SLR = single_stored_object[jj][:SLRModel_SLR])

        df = vcat(df, df_to_add)
    end

    return df
end

# RCP 4.5 w/ TPs
df_SLR_tipall_RCP45 = extract_global_SLR(string(dir_meta_results, "/jld2_files/500_RCP45SSP2_persistDistribution_tipall_tdamageBHM.jld2"))
CSV.write(string(dir_meta_results, "/SLR/SLR_n500_RCP45SSP2_tipall.csv"), df_SLR_tipall_RCP45)

# RCP 4.5 w/o TPs
df_T_AT_tipnone_RCP45 = extract_global_SLR(string(dir_meta_results, "/jld2_files/500_RCP45SSP2_persistDistribution_tipnone_tdamageBHM.jld2"))
CSV.write(string(dir_meta_results, "/SLR/SLR_n500_RCP45SSP2_tipnone.csv"), df_T_AT_tipnone_RCP45)

# RCP 8.5 w/ TPs
df_SLR_tipall_RCP85 = extract_global_SLR(string(dir_meta_results, "/jld2_files/500_RCP85SSP2_persistDistribution_tipall_tdamageBHM.jld2"))
CSV.write(string(dir_meta_results, "/SLR/SLR_n500_RCP85SSP2_tipall.csv"), df_SLR_tipall_RCP85)

# RCP 8.5 w/o TPs
df_T_AT_tipnone_RCP85 = extract_global_SLR(string(dir_meta_results, "/jld2_files/500_RCP85SSP2_persistDistribution_tipnone_tdamageBHM.jld2"))
CSV.write(string(dir_meta_results, "/SLR/SLR_n500_RCP85SSP2_tipnone.csv"), df_T_AT_tipnone_RCP85)

# RCP 2.6 w/ TPs
df_SLR_tipall_RCP26 = extract_global_SLR(string(dir_meta_results, "/jld2_files/500_RCP3-PD26SSP2_persistDistribution_tipall_tdamageBHM.jld2"))
CSV.write(string(dir_meta_results, "/SLR/SLR_n500_RCP3-PD26SSP2_tipall.csv"), df_SLR_tipall_RCP26)

# RCP 2.6 w/o TPs
df_T_AT_tipnone_RCP26 = extract_global_SLR(string(dir_meta_results, "/jld2_files/500_RCP3-PD26SSP2_persistDistribution_tipnone_tdamageBHM.jld2"))
CSV.write(string(dir_meta_results, "/SLR/SLR_n500_RCP3-PD26SSP2_tipnone.csv"), df_T_AT_tipnone_RCP26)

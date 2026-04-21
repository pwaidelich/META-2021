using CSV
using DataFrames

dir_output = "Y:/gem/10 Research/10 Projects/GREENFIN/70 Equity VaR Tipping Points/META_results_Aug2024"

subfolders = readdir(dir_output) |> filter(f -> occursin(r"n10000", f) && !occursin(r"n10000\.csv$", f))

for subfolder_jj in subfolders
    
    print(subfolder_jj)

    # Get a list of all CSV files in the directory
    files_slr = readdir(string(dir_output, "/", subfolder_jj, "/SLR")) |> filter(f -> occursin(r"SLR_mc\d+\.csv$", f) )

    # Initialize an empty DataFrame to hold the combined data
    combined_df = DataFrame()

    # Loop through each file, read it, and append it to the combined DataFrame
    for file_jj in files_slr
        df = CSV.read(string(dir_output, "/", subfolder_jj, "/SLR/", file_jj), DataFrame)
        append!(combined_df, df)
    end

    # Save the combined DataFrame as a CSV file
    CSV.write(string(dir_output, "/", subfolder_jj, "/SLR.csv"), combined_df)
end

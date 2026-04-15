rm(list = ls())

library(tidyverse)
library(vroom)

dir_output <- "Y:/gem/10 Research/10 Projects/GREENFIN/70 Equity VaR Tipping Points/META_results_Aug2024"

mc_samplesize <- 10000

dir_to_compile <- list.files(dir_output) %>% str_subset(paste0("^n", mc_samplesize))

start <- Sys.time()
dir_selected <- dir_to_compile[2]
variable_selected <- "lossfactor"

if(!variable_selected %in% c("T_AT", "SLR", "lossfactor")) stop("Issue")

paths_selected <- file.path(dir_output, dir_selected,
                       ifelse(variable_selected == "lossfactor", "lossfactors_conspc", variable_selected),
                       paste0("variable_selected_mc", 1:mc_samplesize, ".csv"))

if(length(paths_selected) != mc_samplesize) stop("Unexpected # of files")

if(mean(map_lgl(paths_selected, file.exists)) != 1) stop("Not all MC files exist")


df_out <- vroom::vroom(paths_selected, show_col_types = F, col_types = "d")
write_csv(df_out, file.path(dir_output, dir_selected, paste0(variable_selected, ".csv")))
end <- Sys.time()

end - start
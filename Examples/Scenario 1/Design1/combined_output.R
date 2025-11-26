library(dplyr)

files <- list.files(pattern = "assurance_batch_.*\\.rds$")
batch_list <- lapply(files, readRDS)

# If each batch returns a data frame:
combined <- bind_rows(batch_list)

# Save a single combined file:
saveRDS(combined, "combined_assurance_results.rds")

# Summary assurance:
overall_assurance <- mean(combined$Final_Decision == "Successful")
overall_assurance

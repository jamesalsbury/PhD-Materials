library(jsonlite)
library(dplyr)
library(purrr)

# find all json files matching your pattern
files <- list.files(pattern = "assurance_batch_.*\\.json$")

# read them all
results <- map(files, read_json)

# extract assurance (each is a vector of length 1)
batch_assurance <- map_dbl(results, ~ .x$assurance[[1]])

# compute final combined assurance
final_assurance <- mean(batch_assurance)

cat("Combined Assurance Estimate:\n")
print(final_assurance)

# If you want a bootstrap-like CI across batches:
combined_CI <- quantile(batch_assurance, probs = c(0.025, 0.975))
cat("\nEmpirical 95% CI across batches:\n")
print(combined_CI)

# ============================================================
# WEEK 4 ANALYSIS
# ============================================================

# This script combines the results from the 10 subsamples.
#
# Question 1:
# Calculate the final mean count and statistical uncertainty
# for each bacterial strain.
#
# Question 2:
# Calculate WT-mutant differences separately in each subsample,
# then calculate the mean difference, its statistical uncertainty,
# and its significance.


# ============================================================
# SELECT THE 10 RESULT FILES
# ============================================================

cat("Please select the 10 result CSV files.\n")

file_paths <- choose.files(
  multi = TRUE
)

# Check that exactly 10 files were selected.
if (length(file_paths) != 10) {
  stop("Please select exactly 10 result files.")
}


# ============================================================
# BACTERIAL IDs
# ============================================================

allowed_IDs <- c(
  211, -211,
  321, -321,
  2212, -2212,
  3122, -3122,
  3312, -3312,
  3334, -3334
)

positive_IDs <- c(
  211,
  321,
  2212,
  3122,
  3312,
  3334
)

negative_IDs <- c(
  -211,
  -321,
  -2212,
  -3122,
  -3312,
  -3334
)


# ============================================================
# READ THE 10 SUBSAMPLES
# ============================================================

# Matrix that will store the mean count for every bacterial
# strain in every subsample.
subsample_means <- matrix(
  0,
  nrow = length(file_paths),
  ncol = length(allowed_IDs)
)

colnames(
  subsample_means
) <- allowed_IDs


# Store the number of non-empty events in each subsample.
event_counts <- numeric(
  length(file_paths)
)


for (i in 1:length(file_paths)) {

  results <- read.csv(
    file_paths[i]
  )

  # Check that all expected bacterial IDs are present.
  if (!all(
    allowed_IDs %in% results$ID
  )) {
    stop(
      "One of the files does not contain all expected bacterial IDs."
    )
  }

  # Check that the event_count column is present.
  if (!"event_count" %in% names(results)) {
    stop(
      "One of the files does not contain an event_count column."
    )
  }

  # Make sure bacterial IDs are always in the same order.
  results <- results[
    match(
      allowed_IDs,
      results$ID
    ),
  ]

  # Store the mean count from this subsample.
  subsample_means[i, ] <- results$mean_count

  # Store the number of non-empty events.
  event_counts[i] <- results$event_count[1]
}


# ============================================================
# QUESTION 1:
# FINAL MEAN COUNT FOR EACH BACTERIAL STRAIN
# ============================================================

# Calculate the weighted average across the 10 subsamples.
# The weights are the numbers of non-empty events in each file.

final_means <- numeric(
  length(allowed_IDs)
)

for (i in 1:length(allowed_IDs)) {

  final_means[i] <- weighted.mean(
    subsample_means[, i],
    w = event_counts
  )
}

names(
  final_means
) <- allowed_IDs


# Statistical uncertainty is the standard deviation
# of the 10 subsample means.

final_uncertainties <- apply(
  subsample_means,
  2,
  sd
)


# Create final results table for Question 1.

final_results <- data.frame(
  ID = allowed_IDs,
  mean = final_means,
  uncertainty = final_uncertainties
)


cat(
  "\nFinal average bacterial counts:\n\n"
)

print(
  final_results
)


# ============================================================
# QUESTION 2:
# WT-MUTANT ASYMMETRY
# ============================================================

cat(
  "\nAsymmetry analysis:\n\n"
)


# Empty table to store final asymmetry results.
asymmetry_results <- data.frame()


for (i in 1:length(positive_IDs)) {

  positive_ID <- positive_IDs[i]
  negative_ID <- negative_IDs[i]


  # Find the columns corresponding to the WT and mutant.
  positive_column <- which(
    allowed_IDs == positive_ID
  )

  negative_column <- which(
    allowed_IDs == negative_ID
  )


  # ----------------------------------------------------------
  # DIFFERENCE-FIRST METHOD
  # ----------------------------------------------------------

  # Calculate WT-mutant separately in every subsample.
  #
  # This produces 10 differences:
  # D1, D2, ..., D10

  per_file_difference <- (
    subsample_means[, positive_column] -
    subsample_means[, negative_column]
  )


  # Final asymmetry:
  # mean of the 10 WT-mutant differences.

  mean_difference <- mean(
    per_file_difference
  )


  # Statistical uncertainty:
  # standard deviation of the 10 WT-mutant differences.

  difference_uncertainty <- sd(
    per_file_difference
  )


  # Significance:
  # how many standard deviations the mean difference
  # lies away from zero.

  significance <- abs(
    mean_difference
  ) / difference_uncertainty


  # Save the result.

  asymmetry_results <- rbind(
    asymmetry_results,
    data.frame(
      positive_ID = positive_ID,
      negative_ID = negative_ID,
      mean_difference = mean_difference,
      uncertainty = difference_uncertainty,
      significance = significance
    )
  )


  # Print the result.

  cat(
    "IDs ",
    positive_ID,
    " vs ",
    negative_ID,
    "\n",
    sep = ""
  )

  cat(
    "  Mean WT-mutant difference: ",
    round(
      mean_difference,
      6
    ),
    "\n",
    sep = ""
  )

  cat(
    "  Statistical uncertainty: ",
    round(
      difference_uncertainty,
      6
    ),
    "\n",
    sep = ""
  )

  cat(
    "  Significance: ",
    round(
      significance,
      2
    ),
    " sigma\n",
    sep = ""
  )


  if (significance > 3) {

    cat(
      "  Significant asymmetry (> 3 sigma).\n\n"
    )

  } else {

    cat(
      "  No significant asymmetry (< 3 sigma).\n\n"
    )
  }
}


# ============================================================
# SAVE FINAL RESULTS
# ============================================================

write.csv(
  final_results,
  "final_average_results.csv",
  row.names = FALSE
)

write.csv(
  asymmetry_results,
  "final_asymmetry_results.csv",
  row.names = FALSE
)


cat(
  "\nFinal average results saved as ",
  "'final_average_results.csv'.\n",
  sep = ""
)

cat(
  "Final asymmetry results saved as ",
  "'final_asymmetry_results.csv'.\n",
  sep = ""
)


# ============================================================
# FINISHED
# ============================================================

cat(
  "\nAnalysis complete.\n"
)
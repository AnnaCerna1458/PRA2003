# ============================================================
# ANALYSE ONE BACTERIAL DATA SET
# ============================================================

# This script analyses one raw bacterial data file at a time. 
#
# It:
# 1. Reads all events in the file.
# 2. Selects the 12 bacterial strains used in the project.
# 3. Calculates the momentum magnitude of each selected bacterium.
# 4. Counts the number of bacteria of each strain.
# 5. Calculates the average count per non-empty event.
# 6. Divides bacteria into five momentum bins.
# 7. Calculates WT-mutant differences within each momentum bin.
# 8. Saves the overall and momentum-dependent results as CSV files.
#
# The output from 10 independently analysed data sets is later
# combined using the subsampling method.


# ============================================================
# SELECT INPUT FILE
# ============================================================

cat(
  "Please select the input file (for example output-Set1.txt).\n"
)

file_path <- file.choose()


# Read the selected file.
data <- readLines(
  file_path
)


# Check that the selected file is not empty.
if (length(data) == 0) {

  stop(
    "The selected file is empty."
  )
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


# WT / positive IDs.
positive_IDs <- c(
  211,
  321,
  2212,
  3122,
  3312,
  3334
)


# Corresponding mutant / negative IDs.
negative_IDs <- c(
  -211,
  -321,
  -2212,
  -3122,
  -3312,
  -3334
)


# ============================================================
# MOMENTUM BINS
# ============================================================

# Five momentum bins used in the final analysis:
#
# 1: 0-0.654 GeV/c
# 2: 0.654-1.32 GeV/c
# 3: 1.32-2.68 GeV/c
# 4: 2.68-5.88 GeV/c
# 5: >5.88 GeV/c

momentum_breaks <- c(
  0,
  0.654,
  1.32,
  2.68,
  5.88,
  Inf
)


bin_labels <- c(
  "0-0.654",
  "0.654-1.32",
  "1.32-2.68",
  "2.68-5.88",
  ">5.88"
)


number_of_bins <- 5


# ============================================================
# COUNTERS
# ============================================================

# Counter for the total number of bacteria of each strain.
ID_counts <- numeric(
  length(allowed_IDs)
)

names(
  ID_counts
) <- allowed_IDs


# Matrix for bacterial counts in each momentum bin.
#
# Rows = bacterial IDs
# Columns = momentum bins

momentum_counts <- matrix(
  0,
  nrow = length(allowed_IDs),
  ncol = number_of_bins
)


rownames(
  momentum_counts
) <- allowed_IDs


colnames(
  momentum_counts
) <- bin_labels


# ============================================================
# MOMENTUM FUNCTION
# ============================================================

calculate_momentum <- function(
  px,
  py,
  pz
) {

  momentum_squared <- (
    px^2 +
    py^2 +
    pz^2
  )


  # Momentum squared should never be negative.
  if (momentum_squared < 0) {

    stop(
      "The value under the square root cannot be negative."
    )
  }


  momentum_magnitude <- sqrt(
    momentum_squared
  )


  return(
    momentum_magnitude
  )
}


# ============================================================
# EVENT COUNTER
# ============================================================

# Only non-empty events are used in the denominator.
event_count <- 0


# Start reading at the first line.
line_index <- 1


# ============================================================
# LOOP THROUGH ALL EVENTS
# ============================================================

while (
  line_index <= length(data)
) {

  # ----------------------------------------------------------
  # READ EVENT HEADER
  # ----------------------------------------------------------

  header <- data[
    line_index
  ]


  # Split the header into its two values.
  header_split <- strsplit(
    trimws(header),
    "\\s+"
  )[[1]]


  header_numeric <- as.numeric(
    header_split
  )


  # Each event header should contain:
  # event ID and number of bacteria.
  if (
    length(header_numeric) != 2
  ) {

    stop(
      "Each event header should contain an event ID and the number of bacteria."
    )
  }


  event_ID <- header_numeric[1]

  bacteria_number <- header_numeric[2]


  # ----------------------------------------------------------
  # FIND BACTERIAL ROWS BELONGING TO THIS EVENT
  # ----------------------------------------------------------

  start_bacteria <- (
    line_index + 1
  )


  end_bacteria <- (
    line_index +
    bacteria_number
  )


  # Check that all expected bacterial rows exist.
  if (
    end_bacteria > length(data)
  ) {

    stop(
      "The event header specifies more bacterial rows than are available in the file."
    )
  }


  # ----------------------------------------------------------
  # PROCESS NON-EMPTY EVENTS
  # ----------------------------------------------------------

  if (
    bacteria_number > 0
  ) {

    # Count this as a non-empty event.
    event_count <- (
      event_count + 1
    )


    event_bacteria <- data[
      start_bacteria:end_bacteria
    ]


    # Loop through every bacterium in this event.
    for (
      i in 1:length(event_bacteria)
    ) {

      # Split the row into:
      # px, py, pz and bacterial ID.
      bacterium <- strsplit(
        trimws(
          event_bacteria[i]
        ),
        "\\s+"
      )[[1]]


      bacterium_numeric <- as.numeric(
        bacterium
      )


      # Every bacterial row should have four values.
      if (
        length(bacterium_numeric) != 4
      ) {

        stop(
          "Bacterial rows should contain px, py, pz and bacterial ID."
        )
      }


      # Assign values.
      px <- bacterium_numeric[1]

      py <- bacterium_numeric[2]

      pz <- bacterium_numeric[3]

      code <- bacterium_numeric[4]


      # ------------------------------------------------------
      # ONLY ANALYSE BACTERIA USED IN THIS PROJECT
      # ------------------------------------------------------

      if (
        code %in% allowed_IDs
      ) {

        # Increase overall counter.
        ID_counts[
          as.character(code)
        ] <- (
          ID_counts[
            as.character(code)
          ] + 1
        )


        # Calculate momentum magnitude.
        momentum <- calculate_momentum(
          px,
          py,
          pz
        )


        # Determine which of the five momentum bins
        # this bacterium belongs to.
        bin_index <- findInterval(
          momentum,
          momentum_breaks
        )


        # Check that the bin assignment is valid.
        if (
          bin_index < 1 ||
          bin_index > number_of_bins
        ) {

          stop(
            "A momentum value was not assigned to one of the five bins."
          )
        }


        # Increase the corresponding ID/bin counter.
        momentum_counts[
          as.character(code),
          bin_index
        ] <- (
          momentum_counts[
            as.character(code),
            bin_index
          ] + 1
        )
      }
    }
  }


  # ----------------------------------------------------------
  # MOVE TO NEXT EVENT
  # ----------------------------------------------------------

  line_index <- (
    end_bacteria + 1
  )
}


# ============================================================
# PART 1: OVERALL RESULTS
# ============================================================

# Average number of bacteria per non-empty event.
mean_counts <- (
  ID_counts /
  event_count
)


# Per-file Poisson uncertainty.
#
# This is useful as a diagnostic for one data set.
#
# The FINAL statistical uncertainty used in the project is
# obtained later from the spread of the 10 subsamples,
# not from this Poisson uncertainty.

poisson_uncertainties <- (
  sqrt(ID_counts) /
  event_count
)


overall_results <- data.frame(

  ID = allowed_IDs,

  total_count = as.numeric(
    ID_counts
  ),

  event_count = event_count,

  mean_count = as.numeric(
    mean_counts
  ),

  poisson_uncertainty = as.numeric(
    poisson_uncertainties
  )
)


# ============================================================
# SHOW OVERALL RESULTS
# ============================================================

cat(
  "\nNumber of non-empty events:",
  event_count,
  "\n"
)


cat(
  "\nOverall results:\n\n"
)


print(
  overall_results
)


# ============================================================
# PART 2: MOMENTUM-DEPENDENT RESULTS
# ============================================================

# Convert bacterial counts in each momentum bin into
# average counts per non-empty event.

momentum_means <- (
  momentum_counts /
  event_count
)


# ============================================================
# CALCULATE WT-MUTANT DIFFERENCE IN EACH BIN
# ============================================================

pair_results <- data.frame()


for (
  i in 1:length(positive_IDs)
) {

  positive_ID <- positive_IDs[i]

  negative_ID <- negative_IDs[i]


  for (
    bin_index in 1:number_of_bins
  ) {

    positive_count <- momentum_counts[
      as.character(positive_ID),
      bin_index
    ]


    negative_count <- momentum_counts[
      as.character(negative_ID),
      bin_index
    ]


    positive_mean <- momentum_means[
      as.character(positive_ID),
      bin_index
    ]


    negative_mean <- momentum_means[
      as.character(negative_ID),
      bin_index
    ]


    # Difference-first method:
    #
    # WT-mutant difference is calculated within this
    # individual data set and momentum bin.
    #
    # The 10 file-level differences are later combined
    # using the subsampling method.

    difference <- (
      positive_mean -
      negative_mean
    )


    pair_results <- rbind(
      pair_results,

      data.frame(

        positive_ID = positive_ID,

        negative_ID = negative_ID,

        bin_index = bin_index,

        bin_low = momentum_breaks[
          bin_index
        ],

        bin_high = momentum_breaks[
          bin_index + 1
        ],

        bin_label = bin_labels[
          bin_index
        ],

        positive_count = positive_count,

        negative_count = negative_count,

        positive_mean = positive_mean,

        negative_mean = negative_mean,

        difference = difference
      )
    )
  }
}


# ============================================================
# SHOW MOMENTUM RESULTS
# ============================================================

cat(
  "\nMomentum-dependent WT-mutant differences:\n\n"
)


print(
  pair_results
)


# ============================================================
# SAVE RESULTS
# ============================================================

# Use the name of the original input file in the output names.

input_name <- tools::file_path_sans_ext(
  basename(
    file_path
  )
)


# Overall results.
overall_output_path <- file.path(

  dirname(
    file_path
  ),

  paste0(
    input_name,
    "_results.csv"
  )
)


write.csv(

  overall_results,

  overall_output_path,

  row.names = FALSE
)


# Momentum-dependent results.
momentum_output_path <- file.path(

  dirname(
    file_path
  ),

  paste0(
    input_name,
    "_momentum_results.csv"
  )
)


write.csv(

  pair_results,

  momentum_output_path,

  row.names = FALSE
)


# ============================================================
# FINISHED
# ============================================================

cat(
  "\nOverall results saved to:\n",
  overall_output_path,
  "\n"
)


cat(
  "\nMomentum results saved to:\n",
  momentum_output_path,
  "\n"
)


cat(
  "\nAnalysis complete.\n"
)
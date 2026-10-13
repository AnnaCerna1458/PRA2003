import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


# ------------------------------------------------------------
# FILES TO ANALYSE
# ------------------------------------------------------------

file_paths = [
    "output-Set1_momentum_results.csv",
    "output-Set2_momentum_results.csv",
    "output-Set3_momentum_results.csv",
    "output-Set4_momentum_results.csv",
    "output-Set5_momentum_results.csv",
    "output-Set6_momentum_results.csv",
    "output-Set7_momentum_results.csv",
    "output-Set8_momentum_results.csv",
    "output-Set9_momentum_results.csv",
    "output-Set10_momentum_results.csv"
]


# ------------------------------------------------------------
# BACTERIAL PAIR NAMES
# ------------------------------------------------------------

pair_names = {
    (211, -211): "E. coli",
    (321, -321): "Bacillus subtilis",
    (2212, -2212): "Pseudomonas aeruginosa",
    (3122, -3122): "Streptococcus pneumoniae",
    (3312, -3312): "Mycobacterium tuberculosis",
    (3334, -3334): "Salmonella enterica"
}


# ------------------------------------------------------------
# READ ALL 10 SUBSAMPLES
# ------------------------------------------------------------

all_data = []

for file_number in range(len(file_paths)):

    data = pd.read_csv(
        file_paths[file_number]
    )

    # Identify which subsample the data came from.
    data["subsample"] = file_number + 1

    all_data.append(data)


# Combine all 10 subsamples into one table.
combined_data = pd.concat(
    all_data,
    ignore_index=True
)


# ------------------------------------------------------------
# CALCULATE FINAL ASYMMETRY IN EACH MOMENTUM BIN
# ------------------------------------------------------------

# Each R output already contains:
#
# difference = WT mean - mutant mean
#
# for one file and one momentum bin.
#
# Here we combine the 10 file-level differences.
#
# Final asymmetry = mean of the 10 differences
# Statistical uncertainty = SD of the 10 differences

final_results = (
    combined_data
    .groupby(
        [
            "positive_ID",
            "negative_ID",
            "bin_index",
            "bin_low",
            "bin_high",
            "bin_label"
        ]
    )
    .agg(
        mean_difference=(
            "difference",
            "mean"
        ),
        uncertainty=(
            "difference",
            "std"
        )
    )
    .reset_index()
)


# ------------------------------------------------------------
# CALCULATE SIGNIFICANCE WITHIN EACH MOMENTUM BIN
# ------------------------------------------------------------

# Significance =
# absolute mean WT-mutant difference / SD of the differences
#
# This uses the same definition as the overall asymmetry
# calculation from Week 4.

final_results["significance"] = (
    abs(
        final_results["mean_difference"]
    )
    /
    final_results["uncertainty"]
)


final_results["interpretation"] = (
    final_results["significance"]
    .apply(
        lambda x:
        "significant"
        if x > 3
        else "compatible with zero"
    )
)


# ------------------------------------------------------------
# SAVE FINAL BIN RESULTS
# ------------------------------------------------------------

final_results.to_csv(
    "final_momentum_results.csv",
    index=False
)


print(
    "\nFinal momentum-bin results:\n"
)

print(
    final_results[
        [
            "positive_ID",
            "negative_ID",
            "bin_label",
            "mean_difference",
            "uncertainty",
            "significance",
            "interpretation"
        ]
    ]
)


# ------------------------------------------------------------
# TEST WHETHER ASYMMETRY CHANGES ACROSS MOMENTUM BINS
# ------------------------------------------------------------

# For every bacterial pair:
#
# 1. Fit a straight line through the five momentum bins
#    separately for each of the 10 subsamples.
#
# 2. This produces 10 slopes.
#
# 3. Calculate the mean slope and SD of the 10 slopes.
#
# 4. Trend significance =
#    absolute mean slope / SD of the 10 slopes.
#
# Because the momentum bins do not have equal physical widths,
# this slope is interpreted as a trend ACROSS BINS,
# not as a physical slope per GeV/c.

trend_rows = []


pairs = (
    combined_data[
        [
            "positive_ID",
            "negative_ID"
        ]
    ]
    .drop_duplicates()
)


for _, pair in pairs.iterrows():

    positive_ID = pair["positive_ID"]
    negative_ID = pair["negative_ID"]

    slopes = []


    for subsample_number in range(1, 11):

        subset = combined_data[
            (
                combined_data["positive_ID"]
                == positive_ID
            )
            &
            (
                combined_data["negative_ID"]
                == negative_ID
            )
            &
            (
                combined_data["subsample"]
                == subsample_number
            )
        ].sort_values(
            "bin_index"
        )


        # Fit:
        #
        # difference =
        # slope * bin_index + intercept

        slope, intercept = np.polyfit(
            subset["bin_index"],
            subset["difference"],
            1
        )

        slopes.append(
            slope
        )


    # Mean trend across the 10 subsamples.
    mean_slope = np.mean(
        slopes
    )


    # Statistical uncertainty from the spread
    # of the 10 slopes.
    slope_uncertainty = np.std(
        slopes,
        ddof=1
    )


    # Significance of the momentum trend.
    trend_significance = (
        abs(mean_slope)
        /
        slope_uncertainty
    )


    trend_rows.append(
        {
            "positive_ID": positive_ID,
            "negative_ID": negative_ID,
            "mean_slope": mean_slope,
            "slope_uncertainty": slope_uncertainty,
            "trend_significance": trend_significance
        }
    )


trend_results = pd.DataFrame(
    trend_rows
)


# ------------------------------------------------------------
# ADD INTERPRETATION OF MOMENTUM TREND
# ------------------------------------------------------------

trend_results["interpretation"] = (
    trend_results["trend_significance"]
    .apply(
        lambda x:
        "significant momentum trend"
        if x > 3
        else "no significant momentum trend"
    )
)


print(
    "\nMomentum-dependence trend results:\n"
)

print(
    trend_results
)


# Save trend results.
trend_results.to_csv(
    "momentum_trend_results.csv",
    index=False
)


# ------------------------------------------------------------
# CREATE FINAL FIGURES
# ------------------------------------------------------------

pairs = (
    final_results[
        [
            "positive_ID",
            "negative_ID"
        ]
    ]
    .drop_duplicates()
)


for _, pair in pairs.iterrows():

    positive_ID = pair["positive_ID"]
    negative_ID = pair["negative_ID"]


    pair_data = final_results[
        (
            final_results["positive_ID"]
            == positive_ID
        )
        &
        (
            final_results["negative_ID"]
            == negative_ID
        )
    ].sort_values(
        "bin_index"
    )


    # Find the trend significance for this pair.
    pair_trend = trend_results[
        (
            trend_results["positive_ID"]
            == positive_ID
        )
        &
        (
            trend_results["negative_ID"]
            == negative_ID
        )
    ].iloc[0]


    trend_significance = (
        pair_trend[
            "trend_significance"
        ]
    )


    # Biological name for the title.
    biological_name = pair_names.get(
        (
            int(positive_ID),
            int(negative_ID)
        ),
        f"{int(positive_ID)} vs {int(negative_ID)}"
    )


    # Five equally spaced positions corresponding
    # to the five momentum bins.
    x_positions = np.arange(
        1,
        len(pair_data) + 1
    )


    # --------------------------------------------------------
    # MAKE FIGURE
    # --------------------------------------------------------

    plt.figure(
        figsize=(8, 5)
    )


    # Mean WT-mutant difference with statistical uncertainty.
    plt.errorbar(
        x_positions,
        pair_data[
            "mean_difference"
        ],
        yerr=pair_data[
            "uncertainty"
        ],
        fmt="o",
        capsize=4
    )


    # Zero line = no WT-mutant asymmetry.
    plt.axhline(
        y=0,
        linestyle="--"
    )


    # Momentum-bin labels.
    plt.xticks(
        x_positions,
        pair_data[
            "bin_label"
        ]
    )


    # Clear axis labels.
    plt.xlabel(
        "Momentum P (GeV/c)"
    )

    plt.ylabel(
        "WT - mutant mean count per event"
    )


    # Biological title.
    plt.title(
        f"{biological_name}: "
        f"WT-mutant difference across momentum bins"
    )


    # --------------------------------------------------------
    # ADD ONE MOMENTUM-TREND RESULT
    # --------------------------------------------------------

    plt.text(
        0.03,
        0.95,
        (
            f"Trend across bins: "
            f"{trend_significance:.1f}σ"
        ),
        transform=plt.gca().transAxes,
        verticalalignment="top"
    )


    plt.tight_layout()


    # --------------------------------------------------------
    # SAVE FIGURE
    # --------------------------------------------------------

    output_name = (
        f"asymmetry_"
        f"{int(positive_ID)}_vs_"
        f"{int(negative_ID)}.png"
    )


    plt.savefig(
        output_name,
        dpi=300
    )


    plt.close()


    print(
        "Saved:",
        output_name
    )


# ------------------------------------------------------------
# PRINT CLEAR SUMMARY
# ------------------------------------------------------------

print(
    "\n"
    "===================================================="
)

print(
    "SUMMARY"
)

print(
    "====================================================\n"
)


for _, row in trend_results.iterrows():

    positive_ID = int(
        row["positive_ID"]
    )

    negative_ID = int(
        row["negative_ID"]
    )

    biological_name = pair_names.get(
        (
            positive_ID,
            negative_ID
        ),
        f"{positive_ID} vs {negative_ID}"
    )


    print(
        biological_name
    )

    print(
        "  Mean trend across bins:",
        round(
            row["mean_slope"],
            6
        )
    )

    print(
        "  Trend uncertainty:",
        round(
            row["slope_uncertainty"],
            6
        )
    )

    print(
        "  Trend significance:",
        round(
            row["trend_significance"],
            2
        ),
        "sigma"
    )

    print(
        "  Interpretation:",
        row["interpretation"]
    )

    print()


print(
    "All figures and numerical results have been saved."
)
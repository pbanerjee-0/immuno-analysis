gtifrom sklearn.linear_model import LinearRegression
from scipy.interpolate import make_interp_spline
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simpson
from scipy.signal import find_peaks
from numpy import trapezoid
import pandas as pd
import os
import math
import argparse
from openpyxl.styles import Font
from matplotlib.ticker import MaxNLocator
import streamlit as st

# Imaris file is only for one movie, need to account for Imaris export headings
st.title("T-Cell Calcium Spike Threshold Analysis")

filename = st.file_uploader("Upload your Imaris export file", type=["csv", "xlsx"])
if not filename:
    st.info("Upload an exported CSV or Excel file from Imaris to get started.")
    st.stop()

# first row of Imaris export is blank and needs skipped, second row is the one with headers
movie = (
    pd.read_csv(filename, skiprows=1)
    if filename.name.endswith(".csv")
    else pd.read_excel(filename, skiprows=1)
)

# only select the String column names for proper filtering
movie.columns = movie.columns.astype(str)
channels = movie.columns[movie.columns.str.contains("Intensity")].tolist()

if len(channels) < 2:
    st.error(f"Expected at least 2 'Intensity' columns, found: {channels}")
    st.stop()

# index=1 on CFP so the two dropdowns don't start on the same column
yfp_channel = st.selectbox("Select YFP Channel", options=channels, index=0)
cfp_channel = st.selectbox("Select CFP Channel", options=channels, index=1)

if yfp_channel == cfp_channel:
    st.warning("YFP and CFP are set to the same channel.")

movie_id = st.text_input("Enter the Movie ID", value="A")
threshold_multiplier = st.number_input(
    "Threshold multiplier", value=1.4, min_value=0.1, step=0.1
)
qc = st.checkbox("QC: filter top 1% of CFP and YFP intensities before regression")

run = st.button("Run analysis")



def movieRead(reader, movie_name):

    # process contents of the Imaris software export and copy it into a new dataFrame
    movie = reader[(reader["Movie ID"] == movie_name)].copy()

    if movie.empty:
        print(f"Warning: No dynamic parameters or calcium data found for movie {movie_name}'s Imaris export.")

    return movie


# after obtaining spike threshold from regression, find the number of spikes and contacts based on YFP-CFP ratios
def adjust_spike(spike_threshold, movie, cfp_col, yfp_col, apc_col):

    # distance from the threshold lines signifies the distance of a given point from the spike threshold point and aa positive distance indicates a spike
    movie["Distance from Threshold Line"] = movie[yfp_col] - (movie[cfp_col] * spike_threshold)

    # spike value of a timepoint is set to 1 if dist from threshold > 0, else it is 0. spikes are 1's
    movie["Spike"] = (movie["Distance from Threshold Line"] > 0).astype(int)

    # contact value of a timepoint is set to 1 if shortest dist to APC < 1, else it is 0. contacts are 1's
    movie["Contact"] = (movie[apc_col] < 1).astype(int)


# linear regression analysis and plotting of data
def analyze_and_plot_regression(movie_data, movie_name, threshold_multiplier=1.4):
    
    # find the names of the intensity columns. CFP should be first, then YFP, in the Imaris file. if not, analysis will be incorrect
    intensity_col_names = [col for col in movie_data.columns if "Intensity" in col]
    cfp_col_name = cfp_channel
    yfp_col_name = yfp_channel

    # use column names to extract data and convert to np arrays to perform regression
    cfpData = movie_data[cfp_col_name]
    yfpData = movie_data[yfp_col_name]

    x = cfpData.to_numpy().reshape(-1, 1)
    y = yfpData.to_numpy().reshape(-1, 1)

    if len(x) == 0:
        print(f"No data to plot for {movie_name}")
        return

    model = LinearRegression(fit_intercept=False)
    model.fit(x, y)
    baseLine = model.predict(x)

    plt.figure(figsize=(10, 6))  # Set a standard size
    plt.scatter(x, y, color="grey", label="Data Points", alpha=0.1, s=10)
    plt.plot(
        x,
        baseLine,
        color="black",
        label=f"Regression Line (y = {model.coef_[0][0]:.2f}x)",
    )

    # compute and plot spike threshold line
    spikeSlope = model.coef_[0][0]
    spikeLine = threshold_multiplier * (spikeSlope * x)
    plt.plot(
        x,
        spikeLine,
        color="magenta",
        linestyle="--",
        label=f"Spike Threshold (y = {spikeSlope*threshold_multiplier:.2f}x)",
    )

    ax = plt.gca()  # get current axis

    # this ensures about 10 ticks maximum, choosing multiples of 2, 5, or 10 automatically
    ax.xaxis.set_major_locator(MaxNLocator(nbins=10, integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=10, integer=True))

    # dynamic scaling of data
    # add a 5% buffer so data points aren't hugging the edge of the box
    x_max = np.max(x) * 1.05
    y_max = np.max(y) * 1.05
    plt.xlim(0, x_max)
    plt.ylim(0, y_max)

    # add grid to read scalability
    plt.grid(True, linestyle=":", alpha=0.6)

    # stats and formatting
    aboveSpike = (y > spikeLine).sum()
    totalCells = len(y)
    ratio = aboveSpike / totalCells

    plt.xlabel("CFP Intensity")
    plt.ylabel("YFP Intensity")
    plt.title(f"Linear Regression: {movie_name}")
    plt.legend()

    plt.tight_layout()
    plt.show()

    print(f"--- Stats for {movie_name} ---")

    # pass calculated spike threshold on to adjust_spike
    print(
        f"Total Cells: {totalCells} | Spiking: {aboveSpike} | Threshold: {spikeSlope*threshold_multiplier} | Ratio of Spiking Cells to Total Cells: {ratio:.4f}"
    )
    return spikeSlope * threshold_multiplier


# styling
def bold_header(filename):
    from openpyxl import load_workbook

    wb = load_workbook(filename)
    for sheetname in wb.sheetnames:
        ws = wb[sheetname]
        for cell in ws[1]:  # Row 1
            cell.font = Font(bold=True)
    wb.save(filename)


def main():
    parser = argparse.ArgumentParser(
        description="Calculate spike threshold for a specific set of T-Cell data. Accounts for variation in signal intensity for calcium."
    )
    parser.add_argument("filename", help="The raw Excel file exported from Imaris")
    parser.add_argument("--movie_id", default="A", help="ID for the current movie")
    parser.add_argument(
        "--threshold_multiplier",
        type=float,
        default=1.4,
        help="Multiplier for spike threshold (default: 1.4)",
    )
    parser.add_argument("--qc", action="store_true",
                    help="Filter top 1%% of CFP and YFP intensities before regression")


    # command line arguments allow user to choose a file name for a specific mouse in the experiment and a movie ID from that file
    args = parser.parse_args()
    reader_path = args.filename
    movie_id = args.movie_id
    try:
        reader = pd.read_excel(reader_path, skiprows=1)
    except FileNotFoundError:
        print("Error: File not found.")
        return
    if "Movie ID" in reader.columns:
        movie_data = movieRead(reader, movie_id)
    else:
        movie_data = reader.copy()
        movie_data["Movie ID"] = movie_id
    if movie_data.empty:
        return
    intensity_cols = [col for col in movie_data.columns if "Intensity" in col]
    
    if len(intensity_cols) < 2:
        print(
            f"Error: Expected at least 2 'Intensity' columns, found: {intensity_cols}"
        )
        return
    apc_options = [
        col
        for col in movie_data.columns
        if "shortest distance" in col.lower() and "apc" in col.lower()
    ]
    if not apc_options:
        print(
            f"Error: APC column not found. Found these columns instead: {list(movie_data.columns)}"
        )
        return
    apc_col = apc_options[0]

    if args.qc:
        cfp_upper = movie_data[intensity_cols[0]].quantile(0.99)
        yfp_upper = movie_data[intensity_cols[1]].quantile(0.99)
        n_before = len(movie_data)
        movie_data = movie_data[(movie_data[intensity_cols[0]] <= cfp_upper) & (movie_data[intensity_cols[1]] <= yfp_upper)]
        n_after = len(movie_data)
        print(f"QC filter applied: removed {n_before - n_after} outlier rows (top 1% CFP or YFP intensity)")

    final_threshold = analyze_and_plot_regression(movie_data, movie_id, args.threshold_multiplier)

    if final_threshold is None:
        print(f"Error: Regression failed for {movie_id}. No output written.")
        return
    adjust_spike(
        final_threshold, movie_data, intensity_cols[0], intensity_cols[1], apc_col
    )

    output_dir = os.path.dirname(os.path.abspath(args.filename))
    output_path = os.path.join(output_dir, f"{movie_id}_Threshold_Analysis.xlsx")
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        movie_data.to_excel(writer, sheet_name="All", index=False)
        if "Surpass Object" in movie_data.columns:
            for cell_type, group in movie_data.groupby("Surpass Object"):
                group.to_excel(writer, sheet_name=str(cell_type)[:31], index=False)
    bold_header(output_path)
    print(f"Final analysis exported to {output_path}")


# if __name__ == "__main__":
#    main()
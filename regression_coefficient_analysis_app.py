from sklearn.linear_model import LinearRegression
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
from io import BytesIO

# Imaris file is only for one movie, need to account for Imaris export headings
st.title("T-Cell Calcium Spike Threshold Analysis")

filename = st.file_uploader("Upload your Imaris export file")
if not filename:
    st.info("Upload an exported CSV or Excel file from Imaris to get started.")
    st.stop()

# first row of Imaris export is blank and needs skipped, second row is the one with headers
if filename.name.lower().endswith(".csv"):
    movie = pd.read_csv(filename, skiprows=1)
elif filename.name.lower().endswith((".xlsx", ".xls")):
    movie = pd.read_excel(filename, skiprows=1)
else:
    st.error("Please upload a .csv, .xlsx, or .xls file.")
    st.stop()

# only select the String column names for proper filtering
movie.columns = movie.columns.astype(str)
yfp_cfp = movie.columns[movie.columns.str.contains("Intensity")].tolist()
apc_shortest_dist= movie.columns[movie.columns.str.contains("Shortest Distance") & movie.columns.str.contains("APC")].tolist()
track_options = [c for c in movie.columns if "track" in c.lower() and "id" in c.lower()]

if len(yfp_cfp) < 2:
    st.error(f"Expected at least 2 'Intensity' columns, found: {yfp_cfp}")
    st.stop()
if len(apc_shortest_dist) < 1:
    st.error(f"Expected at least 1 'Shortest Distance' column for APC, found: {apc_shortest_dist}")
    st.stop()
if len(track_options) < 1:
    st.error("No Track ID column found in the export.")
    st.stop()

# index=1 on CFP so the two dropdowns don't start on the same column
yfp_channel = st.selectbox("Select YFP Channel", options=yfp_cfp, index=0)
cfp_channel = st.selectbox("Select CFP Channel", options=yfp_cfp, index=1)
track_col = st.selectbox("Select Track ID column", options=track_options, index=0)

if yfp_channel == cfp_channel:
    st.warning("YFP and CFP are set to the same channel.")

movie_id = st.text_input("Enter the Movie ID", value="A")
threshold_multiplier = st.number_input(
    "Threshold multiplier",
    value=1.4,
    min_value=0.1,
    step=0.001,
    format="%.3f",
)

qc = st.checkbox("QC: filter top 1% of CFP and YFP intensities before regression")

run = st.button("Run analysis")



def movieRead(reader, movie_name):

    # process contents of the Imaris software export and copy it into a new dataFrame
    movie = reader[(reader["Movie ID"] == movie_name)].copy()

    if movie.empty:
        print(f"Warning: No dynamic parameters or calcium data found for movie {movie_name}'s Imaris export.")

    return movie


# after obtaining spike threshold from regression, update the number of spikes and contacts based on the new threshold to allow downstream calculation of YFP CFP ratios
def adjust_spike(spike_threshold, movie, cfp_col, yfp_col, apc_col):

    # distance from the threshold lines signifies the distance of a given point from the spike threshold point and aa positive distance indicates a spike
    movie["Distance from Threshold Line"] = movie[yfp_col] - (movie[cfp_col] * spike_threshold)

    # spike value of a timepoint is set to 1 if dist from threshold > 0, else it is 0. spikes are 1's
    movie["Spike"] = (movie["Distance from Threshold Line"] > 0).astype(int)

    # contact value of a timepoint is set to 1 if shortest dist to APC < 1, else it is 0. contacts are 1's
    movie["Contact"] = (movie[apc_col] < 1).astype(int)

    # surfaces that spike while they contact
    movie["Spike and Contact"] = ((movie["Spike"] == 1) & (movie["Contact"] == 1)).astype(int)


# calculate per track ratios for the movie
def calculate_track_ratios(movie, track_col, celltype_col=None):

    # group by cell type as well so the same Track ID in two cell types doesn't get merged
    group_cols = [track_col] + ([celltype_col] if celltype_col else [])

    # For each Track ID, compute the following parameters
    track_stats = movie.groupby(group_cols).agg(
        **{
            "Spiking Cells": ("Spike", "sum"),              # Sum of 1s gives total spiking cells
            "Contacts": ("Contact", "sum"),                 # Sum of 1s gives total contacting cells
            "Spike and Contact": ("Spike and Contact", "sum"),  # Sum of 1s gives total "both" cells
            "Surfaces": ("Spike", "size"),                  # Counting rows gives total surfaces/duration
        }
    ).reset_index()

    # Rename columns
    track_stats = track_stats.rename(columns={track_col: "TrackID"})
    if celltype_col:
        track_stats = track_stats.rename(columns={celltype_col: "Cell Type"})

    # Spike Ratio = Spiking Cells / Total Surfaces
    track_stats["Spike Ratio"] = track_stats["Spiking Cells"] / track_stats["Surfaces"]

    # Contact Ratio = Contacts / Total Surfaces
    track_stats["Contact Ratio"] = track_stats["Contacts"] / track_stats["Surfaces"]

    # Spiking Contact Ratio = (Spike & Contact) / Contacts
    track_stats["Spiking Contact Ratio"] = (
        track_stats["Spike and Contact"]
        .div(track_stats["Contacts"].replace(0, np.nan))
        .fillna(0)  # Handle division by zero
    )

    track_stats["Dataset"] = "BTracks"

    return track_stats

# new method to calculate spike threshold adapted to the lower spread of each individual dataset
def noise_calibrated_threshold_finder(movie_data):

    cfpData = movie_data[cfp_channel]
    yfpData = movie_data[yfp_channel]

    ratio = yfpData/cfpData

    # median of all the YFP-CFP ratios in the movie
    median_ratio=np.median(ratio)
    lower_spread=ratio[ratio<=median_ratio]

    threshold=median_ratio+1.5*lower_spread

 
# old method for finding threshold by fitting YFP to CFP ratios to a regression line for each dataset
def analyze_and_plot_regression(movie_data, movie_name, threshold_multiplier=1.4):
    
    # find the names of the intensity columns. CFP should be first, then YFP, in the Imaris file. if not, analysis will be incorrect

    # use column names to extract data and convert to np arrays to perform regression
    cfpData = movie_data[cfp_channel]
    yfpData = movie_data[yfp_channel]

    x = cfpData.to_numpy().reshape(-1, 1)
    y = yfpData.to_numpy().reshape(-1, 1)

    if len(x) == 0:
        print(f"No data to plot for {movie_name}")
        return

    model = LinearRegression(fit_intercept=False)
    model.fit(x, y)
    baseLine = model.predict(x)

    (fig, ax) = plt.subplots(figsize=(10, 6))  # Set a standard size
    ax.scatter(x, y, color="grey", label="Data Points", alpha=0.1, s=10)
    ax.plot(
        x,
        baseLine,
        color="black",
        label=f"Regression Line (y = {model.coef_[0][0]:.2f}x)",
    )

    # compute and plot spike threshold line
    spikeSlope = model.coef_[0][0]
    spikeLine = threshold_multiplier * (spikeSlope * x)
    ax.plot(
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
    st.pyplot(fig)
    plt.close(fig)

    st.metric(f"--- Stats for {movie_name} ---", value=f"Total Cells: {totalCells}")
    st.metric("Spiking", value=f"{aboveSpike}")
    st.metric("Threshold", value=f"{spikeSlope*threshold_multiplier}")
    st.metric("Ratio", value=f"{ratio:.4f}")

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


if run:
    if "Movie ID" in movie.columns:
        movie_data = movieRead(movie, movie_id)
    else:
        movie_data = movie.copy()
        movie_data["Movie ID"] = movie_id
    if movie_data.empty:
        st.warning(f"No data found for Movie ID {movie_id}.")
        st.stop()

    if qc:
        cfp_upper = movie_data[cfp_channel].quantile(0.99)
        yfp_upper = movie_data[yfp_channel].quantile(0.99)
        n_before = len(movie_data)
        movie_data = movie_data[(movie_data[cfp_channel] <= cfp_upper) & (movie_data[yfp_channel] <= yfp_upper)]
        n_after = len(movie_data)
        st.write(f"QC filter applied: removed {n_before - n_after} outlier rows (top 1% CFP or YFP intensity)")
   

    spike_threshold = analyze_and_plot_regression(movie_data, movie_id, threshold_multiplier)
    if spike_threshold is not None:
        adjust_spike(spike_threshold, movie_data, cfp_channel, yfp_channel, apc_shortest_dist[0])
        st.write("Analysis complete. Here are the first few rows of the processed data:")
        st.dataframe(movie_data.head())
    else:
        st.error("Analysis could not be completed due to insufficient data.")
        st.stop()

    # cell type label for the track ratios, T-Cell Subset if it's in the export, otherwise Surpass Object
    if "T-Cell Subset" in movie_data.columns:
        celltype_col = "T-Cell Subset"
    elif "Surpass Object" in movie_data.columns:
        celltype_col = "Surpass Object"
    else:
        celltype_col = None

    tracks = calculate_track_ratios(movie_data, track_col, celltype_col)
    st.write("Per-track ratios:")
    st.dataframe(tracks)

    buffer = BytesIO()
    with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
        movie_data.to_excel(writer, sheet_name="All", index=False)
        tracks.to_excel(writer, sheet_name="Tracks", index=False)
        if "Surpass Object" in movie_data.columns:
            for cell_type, group in movie_data.groupby("Surpass Object"):
                group.to_excel(writer, sheet_name=str(cell_type)[:31], index=False)

        # bold the header row of every sheet (replaces bold_header)
        for ws in writer.book.worksheets:
            for cell in ws[1]:
                cell.font = Font(bold=True)

    st.download_button(
        label="Download results (Excel)",
        data=buffer.getvalue(),
        file_name=f"{movie_id}_Threshold_Analysis.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    )
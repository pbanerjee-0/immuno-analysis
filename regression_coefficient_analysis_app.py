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

    # measure how far each value in the lower half is from the median ratio
    lower_half=ratio[ratio<=median_ratio]
    deviations = median_ratio - lower_half

    lower_spread = 1.4826 * np.median(deviations)
    lower_spread = max(lower_spread, 0.1)

    threshold = median_ratio + 1.5 * lower_spread

    return threshold, median_ratio




def plot_spikes(movie_data, movie_name, threshold, median_ratio):
    cfpData = movie_data[cfp_channel]
    yfpData = movie_data[yfp_channel]

    valid = np.isfinite(cfpData) & np.isfinite(yfpData)
    x = cfpData[valid].to_numpy()
    y = yfpData[valid].to_numpy()

    if len(x) == 0:
        st.warning(f"No valid data to plot for {movie_name}")
        return

    # same rule as adjust_spike: a surface spikes if YFP is above threshold * CFP
    is_spike = y > threshold * x
    n_spike = int(is_spike.sum())
    n_total = len(y)
    pct_spiking = n_spike / n_total * 100

    # colors (dark pink for spiking points, brighter magenta kept for the threshold line)
    SPIKE_COLOR = "#C2185B"
    LINE_COLOR = "#E91E8C"
    GREY = "#8C8C8C"

    with plt.rc_context({"font.size": 14, "axes.titlesize": 18, "axes.labelsize": 16}):
        fig, ax = plt.subplots(figsize=(10, 7))

        x_max = max(float(np.max(x)) * 1.05, 1.0)
        y_max = max(float(np.max(y)) * 1.05, 1.0)
        x_line = np.linspace(0, x_max, 200)

        # light pink shading over the region above the threshold (the "spiking zone")
        ax.fill_between(x_line, threshold * x_line, y_max, color=SPIKE_COLOR, alpha=0.06, linewidth=0, zorder=0)

        # non-spiking surfaces (grey) first, spiking surfaces (dark pink) on top
        ax.scatter(x[~is_spike], y[~is_spike], s=28, color=GREY, alpha=0.35,
                   edgecolors="none", zorder=2, label=f"Non-spiking (n = {n_total - n_spike})")
        ax.scatter(x[is_spike], y[is_spike], s=36, color=SPIKE_COLOR, alpha=0.85,
                   edgecolors="white", linewidths=0.4, zorder=3, label=f"Spiking (n = {n_spike})")

        # baseline (median ratio) and spike threshold lines
        ax.plot(x_line, median_ratio * x_line, color="black", linewidth=2, zorder=4,
                label=f"Baseline median ratio (y = {median_ratio:.2f}x)")
        ax.plot(x_line, threshold * x_line, color=LINE_COLOR, linestyle="--", linewidth=2.5, zorder=4,
                label=f"Spike threshold (y = {threshold:.2f}x)")

        # percent spiking, clearly visible in the upper left
        ax.text(0.03, 0.97, f"{pct_spiking:.1f}% spiking surfaces\n({n_spike} of {n_total})",
                transform=ax.transAxes, ha="left", va="top", fontsize=18, fontweight="bold",
                color=SPIKE_COLOR, zorder=5,
                bbox=dict(boxstyle="round,pad=0.5", facecolor="white", edgecolor=SPIKE_COLOR, linewidth=2))

        ax.set_xlim(0, x_max)
        ax.set_ylim(0, y_max)
        ax.set_xlabel("CFP intensity")
        ax.set_ylabel("YFP intensity")
        ax.set_title(f"YFP vs CFP intensity: {movie_name}", fontweight="bold")
        ax.grid(True, linestyle=":", alpha=0.5)
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=8))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=8))
        ax.legend(loc="lower right", frameon=True, framealpha=0.95, fontsize=12)

        fig.tight_layout()
        st.pyplot(fig)

        # high-resolution copy for the poster
        png_buffer = BytesIO()
        fig.savefig(png_buffer, format="png", dpi=300, bbox_inches="tight")
        plt.close(fig)

    st.download_button(
        label="Download plot (PNG, 300 dpi)",
        data=png_buffer.getvalue(),
        file_name=f"{movie_name}_spike_plot.png",
        mime="image/png",
    )

    col1, col2, col3 = st.columns(3)
    col1.metric("Spiking Surfaces", f"{n_spike} / {n_total}")
    col2.metric("% Spiking Surfaces", f"{pct_spiking:.2f}%")
    col3.metric("Noise Threshold Ratio", f"{threshold:.3f}")

 
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
   

    # spike_threshold = analyze_and_plot_regression(movie_data, movie_id, threshold_multiplier)
    spike_threshold, median_ratio=noise_calibrated_threshold_finder(movie_data)

    if spike_threshold is not None:
        plot_spikes(movie_data, movie_id, spike_threshold, median_ratio)
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

    # once the spike threshold is locked in, the track ratios can be calculated
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
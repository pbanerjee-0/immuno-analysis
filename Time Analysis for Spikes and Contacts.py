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
from openpyxl.styles import Font
import streamlit as st
from io import BytesIO

# The purpose of this script is to find the duration of every individual calcium spike and contact with APC in the moving using an AUC analysis with scipy.
# The goal is to analyze how frequent calcium spikes are across different datasets, and how the duration of the spike varies.
st.title("T-Cell Calcium Spike and Contact Duration Analysis")

# upload the file exported from the threshold analysis, it needs the Spike and Contact columns already filled in
filename = st.file_uploader("Upload your Threshold Analysis file")
if not filename:
    st.info("Upload the Excel file exported from the threshold analysis to get started.")
    st.stop()

if filename.name.lower().endswith(".csv"):
    reader = pd.read_csv(filename)
elif filename.name.lower().endswith((".xlsx", ".xls")):
    reader = pd.read_excel(filename)
else:
    st.error("Please upload a .csv, .xlsx, or .xls file.")
    st.stop()

reader.columns = reader.columns.astype(str)

for needed in ["Spike", "Contact"]:
    if needed not in reader.columns:
        st.error(f"Expected a '{needed}' column, found: {reader.columns.tolist()}")
        st.stop()

track_options = [c for c in reader.columns if "track" in c.lower() and "id" in c.lower()]
if len(track_options) < 1:
    st.error("No Track ID column found in the file.")
    st.stop()
track_col = st.selectbox("Select Track ID column", options=track_options, index=0)

# cell type is T-Cell Subset if it's in the file, otherwise Surpass Object
celltype_options = [c for c in ["T-Cell Subset", "Surpass Object"] if c in reader.columns]
celltype_col = st.selectbox("Select Cell Type column", options=celltype_options) if celltype_options else None

if "Movie ID" in reader.columns:
    movie_id = st.selectbox("Select the Movie ID", options=reader["Movie ID"].dropna().unique().tolist())
else:
    movie_id = st.text_input("Enter the Movie ID", value="A")

# time between frames in seconds
time = st.number_input("Time between frames (sec)", value=16.466, min_value=0.001, format="%.3f")

run = st.button("Run analysis")


def analyze_durations(movie, track_col, celltype_col, time):

    BTracks = {}
    for trackID, track in movie.groupby(track_col, sort=False):
        spikes   = track['Spike'].tolist()
        contacts = track['Contact'].tolist()
        cellType = track[celltype_col].iloc[0] if celltype_col else ""

        # needed so it doesn't include tracks with 0 duration which leads to float error
        if len(spikes) <= 1:
            continue
        trkDur = time * (len(spikes) - 1)

        spikeDurations=[]
        conDurations=[]

        spikePks=find_peaks(spikes, height=None, width=(None, None))
        numSpikePks=len(spikePks[0])           # find_peaks gives a 2D array, first array is peaks, so this gives # peaks

        # each item in the dictionary that can be retrieved using trackID :
        # # of spiking peaks, # of contacting peaks, track duration, spike frequency, spike duration...and same for contacts
        # last is cell type (CD8 or TReg)

        conPks=find_peaks(contacts, height=None, width=(None, None))
        numConPks=len(conPks[0])
        if len(conPks[0]) == 0:

            # either no contact at all, or it's one big block of 1's
            total_auc = trapezoid(contacts, dx=time)
            conDurations = [total_auc]
            if(conDurations[0]>0):     # if the AUC is more than 1 but length of contact peaks detected is 0, it is one peak
                numConPks=1

        # use left ips and right ips
        # floor the left ips and round up the right ips to get start and end value indices
        for leftips, rightips, w in zip(spikePks[1]['left_ips'], spikePks[1]['right_ips'], spikePks[1]['widths']):

            start=math.floor(leftips)           # spike starts happening from 0 before, so we do this to add in triangle AUC
            end=math.ceil(rightips)

            spikePeakDur=spikes[start : end+1]  # we have to add 1 because the end itself is excluded

            if(len(spikes)==len(spikePks[0])):                   # if the whole track is peaking
                spikeDurations.append(trapezoid(spikes, dx=time))
            else:
                spikeDurations.append(trapezoid(spikePeakDur, dx=time))

        # same calculations as above but for contacts

        for leftips, rightips, w in zip(conPks[1]['left_ips'], conPks[1]['right_ips'], conPks[1]['widths']):

            start=math.floor(leftips)
            end=math.ceil(rightips)

            conPeakDur=contacts[start : end+1]  # we have to add 1 because the end itself is excluded

            if(len(contacts)==sum(contacts)):                   # if the whole track is peaking
                conDurations.append(trapezoid(contacts, dx=time))
            else:
                conDurations.append(trapezoid(conPeakDur, dx=time))

        # plain python floats rounded so they don't show up as np.float64(...) or 402.71999999999997
        spikeDurations = [round(float(d), 4) for d in spikeDurations]
        conDurations   = [round(float(d), 4) for d in conDurations]

        # no coefficients, area under the curve is the duration we want
        BTracks[trackID] = [
            numSpikePks,                    # spiking peaks
            numConPks,                      # contact peaks
            trkDur,                         # track duration
            numSpikePks / (trkDur/3600),    # frequency of spikes, per hour so sec / 3600
            spikeDurations,                 # spike duration
            numConPks   / (trkDur/3600),    # frequency of contacts
            conDurations,                   # contact duration
            cellType,
        ]

    df_B = pd.DataFrame.from_dict(
        BTracks, orient='index',
        columns=[
            'Spiking Peaks','Contact Peaks','Track Duration',
            'Spike Frequency','Spike Durations',
            'Contact Frequency','Contact Durations',
            'Cell Type'
        ]
    ).reset_index().rename(columns={'index':'TrackID'})
    df_B['Dataset'] = 'BTracks'

    return df_B


if run:
    if "Movie ID" in reader.columns:
        movie = reader[reader['Movie ID'] == movie_id].copy()
    else:
        movie = reader.copy()
    if movie.empty:
        st.warning(f"No data found for Movie ID {movie_id}.")
        st.stop()

    df_B = analyze_durations(movie, track_col, celltype_col, time)
    if df_B.empty:
        st.error("No tracks with more than one timepoint were found.")
        st.stop()

    # no lists in the per-track table, the durations get their own sheets with one number per row
    export_B = df_B.drop(columns=['Spike Durations', 'Contact Durations'])

    st.write("Analysis complete. Here are the per-track results:")
    st.dataframe(export_B)

    buffer = BytesIO()
    with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
        export_B.to_excel(writer, sheet_name='Tracks', index=False)

        # every individual duration in one column, one per row, with its track and cell type
        for dur_col in ['Spike Durations', 'Contact Durations']:
            long_df = (
                df_B[['TrackID', 'Cell Type', dur_col]]
                .explode(dur_col)
                .dropna(subset=[dur_col])
                .rename(columns={dur_col: 'Duration'})
            )
            long_df['Duration'] = long_df['Duration'].astype(float)
            long_df.to_excel(writer, sheet_name=dur_col, index=False)

        # bold the header row of every sheet
        for ws in writer.book.worksheets:
            for cell in ws[1]:
                cell.font = Font(bold=True)

    st.download_button(
        label="Download results (Excel)",
        data=buffer.getvalue(),
        file_name=f"{movie_id}_FinalTimeAnalysis.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    )
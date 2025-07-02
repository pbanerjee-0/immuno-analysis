from sklearn.linear_model import LinearRegression
from scipy.interpolate import make_interp_spline
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simpson
from scipy.signal import find_peaks
from numpy import trapz
import pandas as pd
import os
import math


def main():
    reader_path = r'C:\Users\Pami Banerjee\OneDrive - University of Pittsburgh\Desktop\Camirand Lab Research Documents\M796 Matrix.xlsx'
    reader = pd.read_excel(reader_path)

    time = 16.466
    movie = reader[reader['Movie ID'] == 'M796_A']

    BTracks = {}
    for trackID, spike, contact, cellType in zip(
            movie['Track ID'], movie['Spike'], movie['Contact'], movie['T-Cell Subset']
        ):
        if trackID not in BTracks:
            BTracks[trackID] = [0, 0, 0, 0, 0, 0, 0, cellType]
        spikes   = reader[reader['Track ID'] == trackID]['Spike'].tolist()
        contacts = reader[reader['Track ID'] == trackID]['Contact'].tolist()

        spikeDurations=[]
        conDurations=[]

        spikePks=find_peaks(spikes, height=None, width=(None, None))  # figure out how to get these peaks and use indices to get AOC
        numSpikePks=len(spikePks[0])           # find_peaks gives a 2D array, first array is peaks, so this gives # peaks

        # each item in the dictionary that can be retrieved using trackID :
        # # of spiking peaks, # of contacting peaks, track duration, spike frequency, spike duration...and same for contacts
        # last is cell type (CD8 or TReg)

        conPks=find_peaks(contacts, height=None, width=(None, None))
        numConPks=len(conPks[0])
        if(len(spikes)>1):          # needed so it doesn't include tracks with 0 duration which leads to float error
            trkDur   = time * (len(spikes) - 1)
        if(trkDur<1):
            print(trkDur)
        if len(conPks[0]) == 0:

            # either no contact at all, or it's one big block of 1's
            total_auc = trapz(contacts, dx=time)
            conDurations = [total_auc]
            if(conDurations[0]>0):     # if the AUC is more than 1 but length of contact peaks detected is 0, it is one peak
                numConPks=1

        # use left ips and right ips
        # floor the left ips and round up the right ips to get start and end value indices
        for leftips, rightips, w in zip(spikePks[1]['left_ips'], spikePks[1]['right_ips'], spikePks[1]['widths']):
            peak=[]
            width=int(round(w))


            start=math.floor(leftips)           # spike starts happening from 0 before, so we do this to add in triangle AUC
            end=math.ceil(rightips)

            spikePeakDur=spikes[start : end+1]  # we have to add 1 because the end itself is excluded

            if(len(spikes)==len(spikePks[0])):                   # if the whole track is peaking
                spikeDurations.append(trapz(spikes, dx=time))
            else:
                spikeDurations.append(trapz(spikePeakDur, dx=time))

            """
            print(trackID)
            print(spikePeakDur)
            print(width)
            print("Left: ", leftips)
            print("Right: ", rightips)
            print("Start: ", start)
            print("End: ", end)
            """

        # same calculations as above but for contacts

        for leftips, rightips, w in zip(conPks[1]['left_ips'], conPks[1]['right_ips'], conPks[1]['widths']):
            peak=[]
            width=int(round(w))
            """
            print(trackID)
            print(len(contacts))
            print(sum(contacts))
            """

            start=math.floor(leftips)
            end=math.ceil(rightips)

            conPeakDur=contacts[start : end+1]  # we have to add 1 because the end itself is excluded

            if(len(contacts)==sum(contacts)):                   # if the whole track is peaking
                conDurations.append(trapz(contacts, dx=time))
            else:
                conDurations.append(trapz(conPeakDur, dx=time))



        BTracks[trackID][0] = numSpikePks
        BTracks[trackID][1] = numConPks
        BTracks[trackID][2] = trkDur
        BTracks[trackID][3] = numSpikePks / (trkDur/3600)  # frequency of spikes, per hour so sec / 3600
        BTracks[trackID][4] = spikeDurations               # spike duration
        BTracks[trackID][5] = numConPks   / (trkDur/3600)  # frequency of contacts
        BTracks[trackID][6] = conDurations          # change to contact duration

        # no coefficients, area under the curve is the duration we want

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

    output_dir  = os.path.dirname(reader_path)
    output_file = os.path.join(output_dir, 'M796A_FinalTimeAnalysis.xlsx')

    print("Exporting to:", output_file)
    print("DataFrame shape:", df_B.shape)

    try:
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            df_B.to_excel(writer, sheet_name='Tracks', index=False)
        print("✅ Export successful")
    except Exception as e:
        print("❌ Export failed:", e)


if __name__ == '__main__':
    main()


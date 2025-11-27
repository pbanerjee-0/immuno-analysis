from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#Steps to use script for data analysis

#This script can be used for each intravital analysis movie, then the movie data is compiled for the whole mouse in Excel. i.e. - M114

#1. In the Excel spreadsheet that contains the matrix, include the information
    #- T-Cell Subset (E2Crimson POS or NEG)
    # Intensity Median for the CFP and YFP Channels
    # Movie ID (include mouse number and movie letter, preferably)
    # Track ID (Copy and paste from Imaris export, then add a suffix for the movie name and cell-type. ie - "1000000001_ACD8". if this is not done, script will not be accurate
    # Distance from Spike Threshold [ Formula : YFP - (CFP x Slope of Spike Equation) ]
    # Spike (1 if Spike, 0 if No Spike) (Use Distance from Spike Threshold to Determine)
    # Contact (1 if Contact, 0 if No Contact) (Use Shortest Distance from APCs to determine)

#2. Run this script, regression graph will pop up. Take the slope from the Spike Threshold equation.
# This slope will be used in the formula that calculates the distance of each spike from the spike threshold line.
# Use this in the Excel formula, replacing the values using the new threshold slope. Save the sheet with the new values.

#3. Re-run the script. Now open the spreadsheet which is exported with the TrackRatios.
#These per track spike ratios will be pasted into the other spreadsheet.

#4. For the contact analysis, have it done in Excel first. Take the shortest distance to APCs as exported from statistics.
# If it is greater than 1, classify as no contact (0). If distance is less than or equal to 1, classify as contact. (1)


#5. For this script, the only lines you would need to update based on what movie/spreadsheet you use and export are
    # Line 38
    # Line 49
    # Line 127
def main():
    # Read all the data from the Excel file using reader, a pandas object.
    # Paste file path between parentheses to read Excel Sheet
    # This is a line that you would change depending on where you get the Excel data from for the matrix
    reader = pd.read_excel('Input File Name Goes Here')





    #Important : This is one of the only lines you need to change, to read information from only one particular movie in a column
    movie = reader[(reader['Movie ID'] == 'M796_D')] # reader to extract info from any column for the Movie


    #For Spike Analysis -
    # Extract the CFP and YFP columns from the file, write the name of whatever the column headers are
    cfpData = movie['CFP Intensity']
    yfpData = movie['YFP Intensity']
    trackIDs =movie['Track ID']


    # Imports the data from the Spike and Contact columns into lists
    movie_spike=(movie['Spike']).tolist()
    movie_contact=(movie['Contact']).tolist()
    movie_both=(movie['Spike and Contact']).tolist() #Surfaces that spike while they contact


    #Use the movie readers to get the Track ID column, but only for the given movie

    movie_tracks=(movie['Track ID']).tolist()



    BTracks={}


    # Extract cellType for each movie

    movie_cellType = movie['T-Cell Subset'].tolist()

    # Combine trackID, spike, contact, spiking contacts, and cellType into tuples
    allTrackSurfaces= tuple(zip(movie_tracks, movie_spike, movie_contact, movie_both, movie_cellType))


    # Update the dictionaries

    for trackID, spike, contact, both, cellType in allTrackSurfaces:
        if trackID not in BTracks:
            numSpike = 0    # Number of Spiking Surfaces
            numContact=0    # Number of Surfaces that contact APCs
            numBoth=0       # Number of Surfaces that Spike and Contact
            numSurface = 0  # Total Number of Surfaces

            spikeRatio = 0  # Number of Spiking Surfaces / Total Surfaces
            contactRatio=0  # Number of Contact Surfaces / Total Surfaces
            bothRatio=0     # Number of Surfaces that Spike AND Contact / Total Surfaces

            BTracks[trackID] = [numSpike, numContact, numBoth, numSurface, spikeRatio, contactRatio, bothRatio, cellType]  # Add cellType
        BTracks[trackID][0] += spike
        BTracks[trackID][1] += contact
        BTracks[trackID][2] +=both
        BTracks[trackID][3] += 1
        BTracks[trackID][4] = (BTracks[trackID][0]) / (BTracks[trackID][3]) # compute spike ratio for the track
        BTracks[trackID][5] = (BTracks[trackID][1]) / (BTracks[trackID][3]) # compute contact ratio for the track

        # this if-statement is important to avoid zero division error
        if BTracks[trackID][1]!=0:
            BTracks[trackID][6] = (BTracks[trackID][2]) / (BTracks[trackID][1]) # fraction of APC contacts that spike


# Exporting the Ratios to Excel Sheet

        # Function to convert dictionary to DataFrame
    def dict_to_dataframe(track_dict, name):
        df = pd.DataFrame.from_dict(
            track_dict, orient='index',
            columns=['Spiking Cells', 'Contacts', 'Spike and Contact', 'Surfaces', 'Spike Ratio', 'Contact Ratio', 'Spiking Contact Ratio', 'Cell Type']  # Include 'Cell Type'
        )
        df.reset_index(inplace=True)
        df.rename(columns={'index': 'TrackID'}, inplace=True)
        df['Dataset'] = name  # Add an identifier column (optional)
        return df


        #add another column to be exported, which will note down whether the cell is TReg or CD8.
        #need to add T-cell subset as one of the dictionary info bits to ensure this can be exported. need to rewrite code a bit.

    # Convert dictionaries to DataFrames



    df_B = dict_to_dataframe(BTracks, 'BTracks')



    # Export to a single Excel file with multiple sheets
    # Changing this is important, keeps track ratio data in a separate file each time.
    with pd.ExcelWriter('M796DTrackRatios.xlsx', engine='openpyxl') as writer:


        df_B.to_excel(writer, sheet_name='Tracks', index=False)



    #Not necessary to change, but good for clarity
    print("Data exported successfully")



# All Linear Regression and Graphing Stuff Here

    # Convert to numpy arrays so that the linear regression function can be used
    x = cfpData.to_numpy().reshape(-1, 1)   # Set every point from CFP Intensity Column as x
    y = yfpData.to_numpy().reshape(-1, 1)   # Set every point from YFP Intensity Column as y

    # Linear regression - force the line to pass through (0, 0)
    model = LinearRegression(fit_intercept=False)  # Forcing intercept to 0
    model.fit(x, y)

    # Predicting a baseline (no intercept, just slope * x)
    baseLine = model.predict(x)

    # Plot data points and regression line
    plt.scatter(x, y, color='grey', label='Data Points')
    plt.plot(x, baseLine, color='black', label=f'Regression Line (y = {model.coef_[0][0]:.2f}x)')

    # Adding another line to signify spikes
    spikeSlope = model.coef_[0][0]  # Slope
    spikeInter = 0  # Fixed spike threshold intercept
    spikeLine = 1.4*(spikeSlope * x + spikeInter)  # y = slope * x + 1.4
    plt.plot(x, spikeLine, color='magenta', linestyle='--', label=f'Spike Threshold (y = {spikeSlope*1.4:.2f}x + {spikeInter})')


    aboveSpike = (y > spikeLine).sum()

    totalCells = len(y)

    ratio=aboveSpike/totalCells

    # print all counts
    print(f"Total number of cells: {totalCells}")
    print(f"Number of cells above the spike threshold line: {aboveSpike}")
    print(f"Ratio of spiking cells to total cells: {ratio}")

    # Add labels, title, and legend. You can change the title based on the movie you are working on
    plt.xlabel('CFP Intensity')
    plt.ylabel('YFP Intensity')
    plt.title('Linear Regression of CFP vs YFP Data for M796')
    plt.legend()

    # Customizing axis intervals
    plt.xticks(np.arange(0, max(x)[0] + 1, step=2))
    plt.yticks(np.arange(0, max(y)[0] + 1, step=2))

    # Display plot and equation
    plt.ylim(bottom=0)  # Ensure y-axis starts at 0
    plt.show()

    # Print the equation of the forced regression line
    slope = model.coef_[0][0]  # Slope (no intercept)
    print(f"The equation of the regression line in the form y = mx is: y = {slope:.2f}x")




main()

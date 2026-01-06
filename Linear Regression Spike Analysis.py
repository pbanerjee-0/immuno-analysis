from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Load data for a specific movie
def movieRead(reader, movie_name):
    # Filter for the specific movie, using copy to avoid warnings
    movie = reader[(reader['Movie ID'] == movie_name)].copy()

    if movie.empty:
        print(f"Warning: No data found for movie {movie_name}")
        
    return movie

# Calculate per track ratios for each movie
def calculate_track_ratios(movie):
    # For each Track ID, compute the following parameters
    track_stats = movie.groupby('Track ID').agg({
        'Spike': 'sum',               # Sum of 1s gives total spiking cells
        'Contact': 'sum',             # Sum of 1s gives total contacting cells
        'Spike and Contact': 'sum',   # Sum of 1s gives total "both" cells
        'Track ID': 'count',          # Counting rows gives total surfaces/duration
        'T-Cell Subset': 'first'      # Read the label of the T-Cell subset (CD8 or TReg)
    })

    # Rename columns
    track_stats = track_stats.rename(columns={
        'Spike': 'Spiking Cells',
        'Contact': 'Contacts',
        'Spike and Contact': 'Spike and Contact',
        'Track ID': 'Surfaces',
        'T-Cell Subset': 'Cell Type'
    })

    # Calculate all ratios simultaneously
    
    # Spike Ratio = Spiking Cells / Total Surfaces
    track_stats['Spike Ratio'] = track_stats['Spiking Cells'] / track_stats['Surfaces']
    
    # Contact Ratio = Contacts / Total Surfaces
    track_stats['Contact Ratio'] = track_stats['Contacts'] / track_stats['Surfaces']
    
    # Spiking Contact Ratio = (Spike & Contact) / Contacts
    track_stats['Spiking Contact Ratio'] = (
        track_stats['Spike and Contact']
        .div(track_stats['Contacts'])
        .fillna(0)  # Handle division by zero
    )

    # 4. Reset index
    track_stats = track_stats.reset_index()
    
    return track_stats

# Linear Regression Analysis and Plotting of Data
def analyze_and_plot_regression(movie_data, movie_name):
    # Extract data using the passed dataframe
    cfpData = movie_data['CFP Intensity']
    yfpData = movie_data['YFP Intensity']

    # Convert to numpy arrays
    x = cfpData.to_numpy().reshape(-1, 1)
    y = yfpData.to_numpy().reshape(-1, 1)

    # Linear regression - force the line to pass through (0, 0)
    model = LinearRegression(fit_intercept=False)
    model.fit(x, y)

    # Predicting a baseline
    baseLine = model.predict(x)

    # Plot data points and regression line
    plt.scatter(x, y, color='grey', label='Data Points', alpha=0.5) # alpha makes points transparent
    plt.plot(x, baseLine, color='black', label=f'Regression Line (y = {model.coef_[0][0]:.2f}x)')

    # Adding another line to signify spikes
    spikeSlope = model.coef_[0][0]
    spikeInter = 0
    spikeLine = 1.4 * (spikeSlope * x + spikeInter)
    plt.plot(x, spikeLine, color='magenta', linestyle='--', label=f'Spike Threshold (y = {spikeSlope*1.4:.2f}x)')

    # Calculations for print output
    aboveSpike = (y > spikeLine).sum()
    totalCells = len(y)
    ratio = aboveSpike / totalCells

    print(f"--- Regression Stats for {movie_name} ---")
    print(f"Total number of cells: {totalCells}")
    print(f"Number of cells above threshold: {aboveSpike}")
    print(f"Ratio of spiking cells: {ratio:.4f}")

    # Graph Formatting
    plt.xlabel('CFP Intensity')
    plt.ylabel('YFP Intensity')
    plt.title(f'Linear Regression of CFP vs YFP Data for {movie_name}')
    plt.legend()
    
    # Axis handling
    if len(x) > 0:
        plt.xticks(np.arange(0, np.max(x) + 1, step=20)) # Changed step to 20 for readability
        plt.yticks(np.arange(0, np.max(y) + 1, step=20))
    
    plt.ylim(bottom=0)
    plt.show()

    print(f"Equation: y = {spikeSlope:.2f}x")



def main():
    # 1. SETUP
    reader_path = 'Input File Name Goes Here' # <--- UPDATE THIS
    movie_id = 'M796_D'                       # <--- UPDATE THIS
    
    try:
        reader = pd.read_excel(reader_path)
    except FileNotFoundError:
        print("Error: File not found. Check the file path.")
        return

    # 2. LOAD DATA
    movie_data = movieRead(reader, movie_id)
    
    if movie_data.empty:
        return 

    # 3. CALCULATE RATIOS
    df_B = calculate_track_ratios(movie_data)
    df_B['Dataset'] = 'BTracks' 

    # 4. EXPORT TO EXCEL
    output_filename = f'{movie_id}_TrackRatios.xlsx'
    with pd.ExcelWriter(output_filename, engine='openpyxl') as writer:
        df_B.to_excel(writer, sheet_name='Tracks', index=False)
        
    print(f"✅ Data exported successfully to {output_filename}")
    
    # 5. RUN REGRESSION ANALYSIS
    analyze_and_plot_regression(movie_data, movie_id)

if __name__ == "__main__":
    main()
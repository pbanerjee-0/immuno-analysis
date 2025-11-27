import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# The purpose of this script is to create a four quadrant graph that compares T-cells in different conditions : spiking, not spiking, contact with APC, no contact with APC

def main():
    # make a separate excel spreadsheet for each movie, this will allow for a cohesive graph
    reader = pd.read_excel('Input File Name Goes Here')
    contact = reader['Contact']
    spike   = reader['Spike']

    x=contact + np.random.uniform(-0.1, 0.1, size=len(contact))
    y=spike + + np.random.uniform(-0.1, 0.1, size=len(spike))
    

    # 3. Plot
    plt.figure(figsize=(6,6))
    plt.scatter(x, y, color='grey', alpha=0.6)

    # 4. Axis tweaks
    plt.xticks([0,1], ['0','1'])
    plt.yticks([0,1], ['0','1'])
    plt.xlabel('Contact')
    plt.ylabel('Spike')
    plt.title('Contact vs. Spike Scatter Plot for M114')
    plt.grid(True, linestyle='--', alpha=0.3)

    # 5. Show it!
    plt.tight_layout()
    plt.show()


    main()

# Write your code here :-)
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def main():
    # make a separate excel spreadsheet for each movie, this will allow for a cohesive graph
    reader = pd.read_excel('C:\\Users\\Pami Banerjee\\OneDrive - University of Pittsburgh\\Desktop\\Camirand Lab Research Documents\\M114 Cre- Matrix.xlsx')
    contact = reader['Contact']
    spike   = reader['Spike']

    x=contact + np.random.uniform(-0.1, 0.1, size=len(contact))
    y=spike + + np.random.uniform(-0.1, 0.1, size=len(spike))
    """
    Reads an Excel file and makes a scatter plot:
      - X-axis = Contact (0=no, 1=yes)
      - Y-axis = Spike   (0=no, 1=yes)
    """


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

# Example usage:
if __name__ == '__main__':
    path = r'C:\Users\Pami Banerjee\OneDrive - University of Pittsburgh\Desktop\Camirand Lab Research Documents\M795 Matrix Kristina.xlsx'
    main()

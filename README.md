Computational Analysis of Immune Cells that Undergo Various Treatments

During my time as an undergraduate research assistant at the University of Pittsburgh, I have written some scripts to analyze data on T-cells, specifically CD8+, CD4+, and Tregs. 
The goal is to analyze the calcium concentrations of cells. The spikes in calcium and average concentration fluctuates across datasets, which is why the definition of what a calcium spike is
must be accounted for. The linear regression script takes the calcium intensity data, fits the points to a curve, and is used to calculate the spike threshold for each movie : the 

The linear regression based script is used to first find a calcium threshold, then use that calcium threshold as a parameter to calculate the number of calcium spikes in a dataset. The number of calcium spikes 
is essential to find out how extensively immune response is triggered across datasets. The number of spikes and number of contacts with APCs are used to calculate coefficients which are exported
and plotted in GraphPad Prism.

The time analysis script is also used to analyze calcium concentration. However, this uses an area-under-the-curve calculation to derive the duration of each calcium spike in the dataset, as well as the frequency 
of calcium spiking for each track. I have done research on using matplotlib, numpy, pandas, and especially scipy in order to extract the information about the number of calcium peaks and the AUC for
every peak.

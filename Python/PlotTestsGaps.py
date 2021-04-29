import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load in parameters
params = pd.read_csv('/mnt/extraspace/GaiaSelectionFunction/Code/C++/Test/gaps/Optimiser_Properties.dat',skipinitialspace=True)
Nt = params['Nt'][0]

# Load in transformed parameters
xt = pd.read_csv('/mnt/extraspace/GaiaSelectionFunction/Code/C++/Test/gaps/FinalPosition_TransformedParameters.dat',header=None)[0][:Nt].values

# Load in gaps
gaps = pd.read_csv('/mnt/extraspace/GaiaSelectionFunction/TestSets/gaps/edr3_gaps.csv')

# Plot
tbeg, tend = 1717.6256+(np.linspace(1666.4384902198801, 2704.3655735533684, 2) + 2455197.5 - 2457023.5 - 0.25)*4
bins = np.linspace(tbeg,tend,Nt+1)
plt.hist(0.5*(bins[1:]+bins[:-1]),bins = bins,weights=xt,lw=0.3,color='k',histtype='step')
for start,end in zip(gaps['tbeg'].values,gaps['tend'].values):
    plt.fill_between([start,end],y1=[-10,-10],y2=[10,10],color='lightblue')
plt.ylim([-1.05*np.abs(np.max(xt)),+1.05*np.abs(np.max(xt))])
plt.xlabel('OBMT (revolutions)')
plt.ylabel('xt')
plt.savefig('/mnt/extraspace/GaiaSelectionFunction/Code/C++/Test/gaps/plot.png',dpi=300,facecolor='w',bbox_inches='tight')
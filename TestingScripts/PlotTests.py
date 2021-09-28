import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import special
import healpy as hp

# Load in parameters
#directory = '/mnt/extraspace/GaiaSelectionFunction/Code/C++/Test/SlimTest/'
directory = '/Users/douglasboubert/Science/gaia-selection-function/Gaia-EDR3-selection-function/TestingScripts/data/Diagnostic66_activeScaling_quadratic/'
params = pd.read_csv(directory+'OptimiserProperties.dat',skipinitialspace=True,sep=' = ',header=None, index_col=0).T
Nt = int(params['Nt'].values[0])
Nm = int(params['Nm'].values[0])
Nl = int(params['Nl'].values[0])
alpha = 0.5*np.log(2)
a = -3.0

# Load in gaps
gaps = pd.read_csv('/Users/douglasboubert/Science/gaia-selection-function/Gaia-EDR3-selection-function/TestingScripts/data/edr3_gaps.csv')

# Load in transformed parameters
#file = 'TempPosition168_TransformedParameters.dat'
#file = 'TempPositions/TempPosition24_TransformedParameters.dat'
file = 'FinalPosition_TransformedParameters.dat'
xt = pd.read_csv(directory+file,header=None,skipinitialspace=True)[0][:Nt].values
xml = np.reshape(pd.read_csv(directory+file,header=None)[0][Nt:Nt+Nl*Nm].values,(Nl,Nm)).T
pml = np.exp(-2.0*alpha*np.exp(-xml))
pml[xml<a] = np.exp(-2.0*alpha*(1.0-xml[xml<a]+a)*np.exp(-a))

# Plot - xt
plt.figure(figsize=(10,5))
tbeg, tend = 1717.6256+(np.linspace(1666.4384902198801, 2704.3655735533684, 2) + 2455197.5 - 2457023.5 - 0.25)*4
bins = np.linspace(tbeg,tend,Nt+1)
ymin,ymax = -1.05*np.abs(np.max(xt)),+1.05*np.abs(np.max(xt))
plt.hist(0.5*(bins[1:]+bins[:-1]),bins = bins,weights=xt,lw=0.1,color='k',histtype='step')
for start,end in zip(gaps['tbeg'].values,gaps['tend'].values):
    plt.fill_between([start,end],y1=[ymin,ymin],y2=[ymax,ymax],color='lightblue',lw=0.0)
plt.ylim([ymin,ymax])
plt.xlabel('OBMT (revolutions)')
plt.ylabel('xt')
plt.xticks(np.arange(1000,5201,250))
#plt.yticks(np.arange(0.0,1.01,0.1))
plt.grid(axis='both',which='both',lw=0.3,color='lightgrey',zorder=0)
plt.savefig(directory+'plot_xt.png',dpi=300,facecolor='w',bbox_inches='tight')
plt.close()

# Plot - pt
plt.figure(figsize=(10,5))
tbeg, tend = 1717.6256+(np.linspace(1666.4384902198801, 2704.3655735533684, 2) + 2455197.5 - 2457023.5 - 0.25)*4
bins = np.linspace(tbeg,tend,Nt+1)
plt.hist(0.5*(bins[1:]+bins[:-1]),bins = bins,weights=special.expit(xt),lw=0.1,color='k',histtype='step')
for start,end in zip(gaps['tbeg'].values,gaps['tend'].values):
    plt.fill_between([start,end],y1=[0,0],y2=[1,1],color='lightblue',lw=0.0)
plt.ylim([0,1])
plt.xlabel('OBMT (revolutions)')
plt.ylabel('pt')
plt.xticks(np.arange(1000,5251,250))
plt.yticks(np.arange(0.0,1.01,0.1))
plt.grid(axis='both',which='both',lw=0.3,color='lightgrey',zorder=0)
plt.savefig(directory+'plot_pt.png',dpi=300,facecolor='w',bbox_inches='tight')
plt.close()

# Plot - pt - hist
plt.figure(figsize=(10,5))
tbeg, tend = 1717.6256+(np.linspace(1666.4384902198801, 2704.3655735533684, 2) + 2455197.5 - 2457023.5 - 0.25)*4
t = np.linspace(tbeg,tend,Nt+1)
t = 0.5*(t[1:]+t[:-1])
t_bins = np.linspace(tbeg,tend,2693+1)
p_bins = np.linspace(0,1,101)
count = np.histogram2d(t,special.expit(xt),bins=(t_bins,p_bins))[0]
from scipy.ndimage import gaussian_filter
count = gaussian_filter(count,sigma=[0,0])
count /= count.max(axis=1)[:,np.newaxis]
plt.imshow(count.T,extent=([tbeg,tend,0,1]),origin='lower',aspect='auto',cmap='Greys',alpha=1.0)
#plt.plot(0.5*(bins[1:]+bins[:-1]),special.expit(xt),lw=0,marker='o',ms=0.3,mew=0,color='k')
for start,end in zip(gaps['tbeg'].values,gaps['tend'].values):
    plt.fill_between([start,end],y1=[0,0],y2=[1,1],color='lightblue',lw=0.0,alpha=0.5)
plt.ylim([0,1])
plt.xlabel('OBMT (revolutions)')
plt.ylabel('pt')
plt.xticks(np.arange(1000,5251,250))
plt.yticks(np.arange(0.0,1.01,0.1))
plt.grid(axis='both',which='both',lw=0.3,color='lightgrey',zorder=0)
plt.savefig(directory+'plot_pt_hist.png',dpi=300,facecolor='w',bbox_inches='tight')
plt.close()

bins = np.arange(1.7,23.05,0.1)
g = 0.5*(bins[1:]+bins[:-1])[:Nm]
if Nm>1 and Nl>12:
    # Plot - xml
    fig, axes = plt.subplots(ncols=4,nrows=5,figsize=(12,10))
    for i in range(5):
        for j in range(4):
            plt.axes(axes[i,j])
            m = 4*i + j
            try:
                hp.mollview(xml[10*m+2], nest=False, xsize=1000, title=f'g={g[10*m+2]:.2f}', hold=True, min=-10,max=10,cmap='RdBu_r', cbar=True,coord=['C','G'],notext=True)
            except IndexError:
                plt.axis('off')
    plt.savefig(directory+'plot_xml.png',dpi=300,facecolor='w',bbox_inches='tight')
    plt.close()

    # Plot - pml
    fig, axes = plt.subplots(ncols=4,nrows=5,figsize=(12,10))
    for i in range(5):
        for j in range(4):
            plt.axes(axes[i,j])
            m = 4*i + j
            try:
                hp.mollview(pml[10*m+2], nest=False, xsize=1000, title=f'g={g[10*m+2]:.2f}', hold=True, min=0,max=1,cmap='viridis', cbar=True,coord=['C','G'],notext=True)
            except IndexError:
                plt.axis('off')
    plt.savefig(directory+'plot_pml.png',dpi=300,facecolor='w',bbox_inches='tight')
    plt.close()

if False:
    # Plot - xml
    fig, axes = plt.subplots(ncols=3,nrows=2,figsize=(9,4))
    for i in range(2):
        for j in range(3):
            plt.axes(axes[i,j])
            m = 3*i + j
            idx = 180+2*m
            try:
                hp.mollview(xml[idx], nest=False, xsize=2000, title=f'g={g[idx]+0.05:.2f}', hold=True, min=-5,max=5,cmap='RdBu_r', cbar=True,coord=['C','G'],notext=True)
            except IndexError:
                plt.axis('off')
    plt.savefig(directory+'plot_xml_fancy.png',dpi=300,facecolor='w',bbox_inches='tight')
    plt.close()


# Plot - xm
xmin, xmax = -1.05*np.max(np.abs(2*xml)),+1.05*np.max(np.abs(2*xml))
plt.figure(figsize=(6,5))
for l in np.random.choice(range(Nl),100):
    plt.plot(g,2*xml[:,l],color='lightgrey',alpha=0.3,lw=0.5)
plt.ylim([xmin,xmax])
plt.xlabel('G')
plt.ylabel('xm')
plt.xticks(np.arange(2,23.1))
plt.savefig(directory+'plot_xm.png',dpi=300,facecolor='w',bbox_inches='tight')
plt.close()

# Plot - pm
for l in np.random.choice(range(Nl),100):
    plt.plot(g,pml[:,l],color='lightgrey',alpha=0.3,lw=0.5)
plt.ylim([0,1])
plt.xlabel('G')
plt.ylabel('pm')
plt.xticks(np.arange(2,23.1))
plt.savefig(directory+'plot_pm.png',dpi=300,facecolor='w',bbox_inches='tight')
plt.close()

# Plot - xm hist
xmin, xmax = -1.05*np.max(np.abs(2*xml)),+1.05*np.max(np.abs(2*xml))
g_bins = np.arange(1.7,23.05,0.1)
x_bins = np.linspace(xmin,xmax,101)
g = 0.5*(g_bins[1:]+g_bins[:-1])[:Nm]
count = np.histogram2d(np.repeat(g[:,np.newaxis],xml.shape[1],axis=1).flatten(),2.0*xml.flatten(),bins=(g_bins,x_bins))[0]
count /= count.max(axis=1)[:,np.newaxis]
plt.imshow(count.T,extent=([1.7,23.0,xmin,xmax]),origin='lower',aspect='auto',cmap='Greys',alpha=1.0)
plt.ylim([xmin,xmax])
plt.xlabel('G')
plt.ylabel('xm')
plt.xticks(np.arange(2,23.1))
plt.savefig(directory+'plot_xm_hist.png',dpi=300,facecolor='w',bbox_inches='tight')
plt.close()

# Plot - pm hist
g_bins = np.arange(1.7,23.05,0.1)
p_bins = np.linspace(0,1,101)
g = 0.5*(g_bins[1:]+g_bins[:-1])[:Nm]
count = np.histogram2d(np.repeat(g[:,np.newaxis],xml.shape[1],axis=1).flatten(),pml.flatten(),bins=(g_bins,p_bins))[0]
count /= count.max(axis=1)[:,np.newaxis]
plt.imshow(count.T,extent=([1.7,23.0,0,1]),origin='lower',aspect='auto',cmap='Greys',alpha=1.0)
plt.ylim([0,1])
plt.xlabel('G')
plt.ylabel('pm')
plt.xticks(np.arange(2,23.1))
plt.savefig(directory+'plot_pm_hist.png',dpi=300,facecolor='w',bbox_inches='tight')
plt.close()

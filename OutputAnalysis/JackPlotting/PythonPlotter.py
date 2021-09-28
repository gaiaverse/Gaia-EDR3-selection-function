import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import special
import healpy as hp
import os
import imageio
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]})
    
def temporalPlot(pt,Nt,gaps):
    tbeg, tend = 1717.6256+(np.linspace(1666.4384902198801, 2704.3655735533684, 2) + 2455197.5 - 2457023.5 - 0.25)*4
    bins = np.linspace(tbeg,tend,Nt+1)
    centres = 0.5*(bins[1:] + bins[:-1])

    plt.plot(centres,pt)
    
    if type(gaps) != bool:
        for start,end in zip(gaps['tbeg'].values,gaps['tend'].values):
            plt.fill_between([start,end],y1=[-2,-2],y2=[2,2],color='lightblue',lw=0.0)
    plt.ylim([0,1.02])
    plt.xlabel('OBMT (revolutions)')
    plt.ylabel('pt')
    plt.xticks(np.arange(1000,5251,250))
    plt.yticks(np.arange(0.0,1.01,0.1))
    # plt.grid(axis='both',which='both',lw=0.3,color='lightgrey',zorder=0)
    # plt.savefig(directory+'plot_pt.png',dpi=300,facecolor='w',bbox_inches='tight')

def spatialPlot(pml):
    hp.mollview(pml, nest=False, hold=True, min=0,max=1,cmap='viridis', cbar=True,coord=['C','G'],notext=True)

def readFile(directory,file):
    try:
        params = pd.read_csv(directory+'OptimizerProperties.dat',skipinitialspace=True,index_col=0,sep=" =",header=None,engine='python').transpose()
        Nt = int(params['Nt'])
        Nm = int(params['Nm'])
        Nl = int(params['Nl'])
    except:
        params = pd.read_csv(directory+'Optimiser_Properties.dat',skipinitialspace=True,sep=",")
        Nt = int(params['Nt'])
        Nm = int(params['Nm'])
        Nl = int(params['Nl'])
    alpha = 0.5*np.log(2)
    a = -3.0
    if type(file) == int:
        file = "TempPositions/TempPosition" + str(file) + "_TransformedParameters.dat"
    xt = pd.read_csv(directory+file,header=None)[0][:Nt].values
    pt = special.expit(xt)
    xml = np.reshape(pd.read_csv(directory+file,header=None)[0][Nt:Nt+Nm*Nl].values,(Nl,Nm)).T
    pml = np.exp(-2.0*alpha*np.exp(-xml))
    pml[xml<a] = np.exp(-2.0*alpha*(1.0-xml[xml<a]+a)*np.exp(-a))
    return [pt,pml]

def baseComparisonPlot(fig,truept,truepml,pt,pml,gaps):
 
    ##temporal plot
    fig.add_subplot(2,1,1)
    temporalPlot(truept,len(truept),False)
    temporalPlot(pt,len(pt),False)
    
    
    ##spatial plot
    fig.add_subplot(2,2,3)
    spatialPlot(truepml[0])
    fig.add_subplot(2,2,4)
    spatialPlot(pml[0])
   
    
def baseGif(saveName,truedirectory,tempdirectory,start,stop,gap):
    
    # ~ trueFile = "True_TransformedParameters.dat"
    # ~ [truePt,truePml] = readFile(truedirectory,trueFile)
    
    names = []
    # ~ fig = plt.figure(figsize=(8,8))
    for i in range(start,stop,gap):
        print(i)
      
        # ~ [pt,pml] = readFile(directory,i)
        
    
        # ~ baseComparisonPlot(fig,truePt,truePml,pt,pml,False)
        
        tempName = saveName +"_temp" + str(i) + ".png"
        names.append(tempName)
        # ~ plt.savefig(tempName)
        # ~ fig.clf()
    # ~ plt.close(fig)
     
    # ~ saveName = saveName + ".gif"
    # ~ delay = 2
    # ~ with imageio.get_writer(saveName,mode = 'I') as writer:
		
        # ~ for name in names:
            # ~ print(name)
            # ~ image = imageio.imread(name)
            # ~ for k in range(0,delay):
                # ~ writer.append_data(image)
            
    for name in names:
        os.remove(name)
        
dataDir = "../../../Data/TestSets/flat/"
directory = "../../../CodeOutput/TestSets/flat/"
baseGif("test",dataDir,directory,4,101,4)

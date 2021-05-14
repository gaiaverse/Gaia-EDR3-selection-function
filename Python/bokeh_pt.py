import pandas as pd
import numpy as np
from scipy import special
import glob
import os.path

# Load in parameters
#directory = '/mnt/extraspace/GaiaSelectionFunction/Code/C++/Output/SmoothingTest/'
#directory = '/mnt/extraspace/GaiaSelectionFunction/Code/C++/Output/'
#root_directory = './data/'
root_directory = '/mnt/extraspace/GaiaSelectionFunction/Code/C++/Output/'

from bokeh.plotting import figure, output_file, show, curdoc
from bokeh.models import ColumnDataSource, Slider, CustomJS, Button, TextInput
from bokeh.layouts import column, row

tbeg, tend = 1717.6256+(np.linspace(1666.4384902198801, 2704.3655735533684, 2) + 2455197.5 - 2457023.5 - 0.25)*4
data = {'t':[tbeg,tend],'existing':[0.5,0.5]}
source = ColumnDataSource(data=data)
#output_file('pt.html')
p = figure(title='Detection probability with time',plot_width=1200, plot_height=400,
           background_fill_color="#fafafa",tools='reset,xbox_zoom,xpan')
p.y_range.start = 0
p.y_range.end = 1
p.x_range.start = tbeg
p.x_range.end = tend
p.xaxis.axis_label = 'OBMT (revolutions)'
p.yaxis.axis_label = 'Detection probability (pt)'

p.line('t', 'existing', line_width=1,source=source)

# Add a slider to control which year the patches are coloured by
epoch_slider = Slider(start=1, end=2, value=1, step=1, title="Epoch", orientation="horizontal", disabled = True)

# Add the custom javascript to interface the slider with the patches
callback = CustomJS(args=dict(source=source, epoch=epoch_slider), code="""
        source.data['existing'] = source.data[String(epoch.value)]
        source.change.emit();
    """)
    
epoch_slider.js_on_change('value', callback)

def refresh():
    directory = root_directory + text.value
    if os.path.isdir(directory):
        
        params = pd.read_csv(directory+'/Optimiser_Properties.dat',skipinitialspace=True)
        Nt = int(params['Nt'][0])
        
        while os.path.isfile(directory+f'/TempPositions/TempPosition{epoch_slider.end+1}_TransformedParameters.dat'):
            n = epoch_slider.end + 1
            file = directory+f'/TempPositions/TempPosition{n}_TransformedParameters.dat'
            source.data[str(n)] = special.expit(pd.read_csv(file,header=None)[0][:Nt].values)
            epoch_slider.end = n
    else:
        epoch_slider.disabled = True

refresh_button = Button(label='Refresh')
refresh_button.on_click(refresh)

# 
text = TextInput(title="Which folder", value='Enter directory here')
# Set up callbacks
def update_directory(attrname, old, new):
    directory = root_directory + text.value
    if os.path.isdir(directory):

        params = pd.read_csv(directory+'/Optimiser_Properties.dat',skipinitialspace=True)
        Nt = int(params['Nt'][0])

        t = np.linspace(tbeg,tend,Nt+1)
        t = 0.5*(t[1:]+t[:-1])
        data = {'t':t}

        # Load in transformed parameters
        files = glob.glob(directory+'/TempPositions/TempPosition*_TransformedParameters.dat')

        N = 0
        for file in files:
            epoch = file.split('/')[-1].split('_')[0][12:]
            data[epoch] = special.expit(pd.read_csv(file,header=None)[0][:Nt].values)
            N += 1
        data['existing'] = data['1']
        epoch_slider.end = N
        epoch_slider.disabled = False
        
        source.data = data
        #source.change.emit()
    else:
        epoch_slider.disabled = True

text.on_change('value', update_directory)

# Combine plot and slider and output
curdoc().add_root(row(p,column(text,epoch_slider,refresh_button)))
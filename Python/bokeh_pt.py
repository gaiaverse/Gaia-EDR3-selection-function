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
from bokeh.models import ColumnDataSource, Slider, CustomJS, Button, Select, Range1d
from bokeh.layouts import column, row

tbeg, tend = 1717.6256+(np.linspace(1666.4384902198801, 2704.3655735533684, 2) + 2455197.5 - 2457023.5 - 0.25)*4
data = {'t':[tbeg,tend],'existing':[0.5,0.5]}
source = ColumnDataSource(data=data)
#output_file('pt.html')
p = figure(title='Detection probability with time',plot_width=1200, plot_height=400,
           background_fill_color="#fafafa",tools='reset,xbox_zoom,xpan',x_range=Range1d(tbeg,tend,bounds=(tbeg, tend)),y_range=Range1d(0,1,bounds=(0, 1)))
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

def fetch_data():
    directory = root_directory + select.value
    
    params = pd.read_csv(directory+'/Optimiser_Properties.dat',skipinitialspace=True)
    Nt = int(params['Nt'][0])
    
    if epoch_slider.disabled == True:
        n = 0
        data = {'t': np.linspace(tbeg, tend, Nt), 'existing':0.5*np.ones(Nt)}
    else:
        data = dict(source.data)
        n = epoch_slider.end if '2' in data.keys() else 1
        
        
    while os.path.isfile(directory+f'/TempPositions/TempPosition{n+1}_TransformedParameters.dat'):
        file = directory+f'/TempPositions/TempPosition{n+1}_TransformedParameters.dat'
        data[str(n)] = special.expit(pd.read_csv(file,header=None,nrows=Nt)[0].values)
        n += 1
    try:
        data['existing'] = data['1']
        source.data = data
        epoch_slider.end = n
        epoch_slider.disabled = False
        epoch_slider.value = 1
    except KeyError:
        pass
    
fetch_data_button = Button(label='Fetch data')
fetch_data_button.on_click(fetch_data)

def change_directory(attrname, old, new):

    epoch_slider.disabled = True
    directories = [directory.split('/')[-2] for directory in glob.glob(root_directory+'*/')]
    select.options = directories

select = Select(title="Which run to load?", value='', options=['Loading directories...'])
select.on_change('value', change_directory)
change_directory('value',0,1)



# Combine plot and slider and output
curdoc().add_root(row(p,column(select,epoch_slider,fetch_data_button)))
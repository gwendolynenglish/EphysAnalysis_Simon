#Import required packages|
import pickle
import os
import ipywidgets
import csv
from ipywidgets import Layout, HBox, VBox
from IPython.display import display
from load_probe_info import *
from cycleDirectoryLFP import *

##Main path for the data 

inputPath_html = ipywidgets.HTML(value = "<p><b>Path to the data of the experiment:</b><br />Enter the path to the folder (with no '/' at the end) that is hierarchically right above the folders of the recording sessions</p>")
inputPath = ipywidgets.Text(value = "/media/loaloa/gdrive/projects/ephys/data", placeholder = "Enter path for data", disabled = False)
# display(VBox([inputPath_html, inputPath]))

##Main path for the output results and figures 

outputPath_html = ipywidgets.HTML(value = "<p><b>Path for the resulting analysis and figures:</b><br />Enter the path to the folder where all results should be stored </p>")
outputPath = ipywidgets.Text(value = "/media/loaloa/gdrive/projects/ephys/output/LFP_output", placeholder = "Enter path for data", disabled = False)
# display(VBox([outputPath_html, outputPath]))

##Sampling rate
sr = ipywidgets.IntText(value = 32000, disabled = False)
# display(VBox([ipywidgets.HTML(value = "<b> Sampling rate (Hz): </b>"),sr]))

##Probe inf1
pi_html = ipywidgets.HTML(value = "<b>Type of the probe used in the experiment</b>")
pi = ipywidgets.Dropdown(options=['a2x16_10mm_100_500_177', 'a2x16_10mm_50_500_177', 'a1x32_6mm_100_177', 'a4x8_5mm_200_400_177','custom'], 
                   value = 'a1x32_6mm_100_177',  disabled = False)
# display(VBox([pi_html, pi]))

##TimeWindow

tw = ipywidgets.Dropdown(options = [('-50-250', 1), ('-50-500', 2)], disabled = False)
# display(VBox([ipywidgets.HTML(value = "<b>Select the time window for analysis(in ms)</b>"), tw]))

#low_pass_freq

lp = ipywidgets.FloatText(value = 500, disabled = False)
# display(VBox([ipywidgets.HTML(value = "<b> Enter the cutoff frequency of the low pass filter to extract LFP from data (in Hz)"), lp]))


#CSD Analysis

csd = ipywidgets.Checkbox(value = True, disabled = False)
# display(HBox([ipywidgets.HTML(value = "<b> Check if Current Source Density analysis should be completed."), csd]))


p = {} #Parameter dictionary (empty)

#Entering the probe info and electrode mapping into the dictionary
probe_info = load_probe_info(pi.value)
p['shanks'] = probe_info['numShanks']

p['probe_name'] = probe_info['name']
p['nr_of_electrodes'] = probe_info['numTrodes']
p['nr_of_electrodes_per_shank'] = probe_info['numTrodesPerShank']
p['nr_of_shanks'] = p['shanks']
p['bottom_ycoord'] = probe_info['bottom_ycoord']
p['top_ycoord'] = probe_info['top_ycoord']
p['id'] = probe_info['id']
p['sitespacing'] = probe_info['sitespacing']

#Entering the path and file format info into the dictionary
p['inputPath'] = inputPath.value
p['outputPath'] = outputPath.value

#Entering the general parameters into the dictionary
p['sample_rate'] = sr.value
    
#Entering the LFP analysis parameters into the dictionary
if tw.value == 1:
    p['evoked_pre'] = 0.05
    p['evoked_post'] = 0.250
if tw.value == 2:
    p['evoked_pre'] = 0.05
    p['evoked_post'] = 0.5
p['low_pass_freq'] = lp.value
        
p['csd'] = csd.value 
    
if not os.path.exists(outputPath.value + '/AnalysisFiles'):
    os.mkdir(outputPath.value + '/AnalysisFiles')
    
#Saving the dictionary in the pickle file and csv named parametersDict
pickle.dump(p, open((outputPath.value + '/AnalysisFiles/parametersDict.p'), 'wb'))

with open(outputPath.value + '/AnalysisFiles/parametersDict.csv', 'w') as textfile:
    fieldnames = ['Field', 'Value']
    writer = csv.DictWriter(textfile, fieldnames = fieldnames)
    writer.writeheader()
    data = [dict(zip(fieldnames, [k,v])) for k, v in p.items()]
    writer.writerows(data)
    
print(p)

analyzeAllFiles(p)
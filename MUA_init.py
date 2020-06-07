import pickle
import os
import ipywidgets
import csv

from load_probe_info import *

def initialize():
    ### 2) Provide Information for Dictionary 
    ##Main path for the data 
    val = ('<p><b>Path to the data of the experiment:</b><br />Enter the path '
           'to the folder (with no `/` at the end) that is hierarchically right'
           ' above the folders of the recording sessions</p>')
    inputPath_html = ipywidgets.HTML(value = val)
    inp_path = "/run/media/loaloa/Samsung_T5/mGE82838485"
    inputPath = ipywidgets.Text(value = inp_path, 
                                placeholder = "Enter path for data", 
                                disabled = False)

    ##Main path for the output results and figures 
    val = ("<p><b>Path for the resulting analysis and figures:</b><br />Enter "
           "the path to the folder (with '/') where all results should be "
           "stored </p>")
    outputPath_html = ipywidgets.HTML(value = val)
    out_path = "/run/media/loaloa/Samsung_T5/output/output"
    outputPath = ipywidgets.Text(value = out_path,
                                placeholder = "Enter path for data", 
                                disabled = False)

    ##Sampling rate
    sr = ipywidgets.IntText(value = 32000, disabled = False)

    ##Probe info
    inp = "<b>Type of the probe used in the experiment</b>"
    pi_html = ipywidgets.HTML(value = inp)
    electrodes = ['a2x16_10mm_100_500_177', 'a2x16_10mm_50_500_177', 
                  'a1x32_6mm_100_177', 'a4x8_5mm_200_400_177']
    pi = ipywidgets.Dropdown(options=electrodes, value = 'a1x32_6mm_100_177',  
                             disabled = False)

    ##TimeWindow
    tw = ipywidgets.Dropdown(options = [('-0.050-0.20', 1), ('-0.010-0.050', 2)], 
                                        value = 1, disabled = False)

    #High_pass_freq
    hp = ipywidgets.FloatText(value = 800, disabled = False)

    #Thresold

    th = ipywidgets.FloatText(value = 8, disabled = False)

    #PSTH binsize 
    psth_bs = ipywidgets.Dropdown(options = [('5', 1), ('10',2)], value = 1,  
                                disabled = False)


    # ### 3) Write Dictionary 
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
        
    #Entering the MUA analysis parameters into the dictionary
    if tw.value == 1:
        p['evoked_pre'] = 0.05
        p['evoked_post'] = 0.20
    if tw.value == 2:
        p['evoked_pre'] = 0.01
        p['evoked_post'] = 0.05
    p['high_pass_freq'] = hp.value
    p['threshold'] = th.value 
    if psth_bs.value == 1:
        p['psth_binsize'] = 5 
    if psth_bs.value == 2:
        p['psth_binsize'] = 10 
            
    if not os.path.exists(outputPath.value + '/AnalysisFiles'):
        os.mkdir(outputPath.value + '/AnalysisFiles')
        
    #Saving the dictionary in the pickle file and csv named parametersDict
    pickle.dump(p, open((outputPath.value+'/AnalysisFiles/parametersDict.p'), 'wb'))

    file = outputPath.value + '/AnalysisFiles/parametersDict.csv'
    with open(file, 'w') as textfile:
        fieldnames = ['Field', 'Value']
        writer = csv.DictWriter(textfile, fieldnames = fieldnames)
        writer.writeheader()
        data = [dict(zip(fieldnames, [k,v])) for k, v in p.items()]
        writer.writerows(data)
    return p
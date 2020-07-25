# coding: utf-8

import numpy as np
import warnings
import concurrent.futures
import os


from  MUA_utility import fetch, slice_data
import MUA_constants as const
from preprocessing import compress_CSVs
import plotting

def process_data(multithreading=7):
    from MUA_cycle_dirs import MUA_analyzeMouseParadigm
    warnings.filterwarnings('ignore')
    
    dirs = os.listdir(const.P['inputPath'])
    # dirs = ['mGE84_30.07.2019_O25C1.mcd']
    if multithreading:
        with concurrent.futures.ProcessPoolExecutor(max_workers=multithreading) as executer:
            [executer.submit(MUA_analyzeMouseParadigm, folder) for folder in dirs]
    else:
        [MUA_analyzeMouseParadigm(folder) for folder in dirs]
    
    warnings.filterwarnings('default')

# process_data()
# compress_CSVs()

# plotting.firingrate_heatmaps('noise_subtr_', subtr_noise='paradigm_wise')
# plotting.firingrate_heatmaps('noisy', False)
# plotting.firingrate_noise_timeline('noisyy')
# plotting.firingrate_noise_timeline(fname_prefix='parad_subtr', subtr_noise='paradigm_wise')

# plotting.make_CSD_summary_plots(lfp_output_appdx='/../LFP_output', dest_dir_appdx='/../CSD_lowthr')
# plotting.plot_time_to_first(dest_dir=f'{const.PROJ_DIR}/output/time_to_first')
# plotting.plot_evoked_lfp(dest_dir_appdx='../thalamus_mapping_lowthr', 
#                          anatomy_dir=f'{const.PROJ_DIR}/metadata/Histology_pngs',
                        #  ts_dir=f'{const.PROJ_DIR}/output/time_to_first')
plotting.oddball10_si(which='O10', dest_dir_appdx='../plots/oddball_si', fname_appdx='100-200ms_Oddball_10')
plotting.oddball10_si(which='O25', dest_dir_appdx='../plots/oddball_si', fname_appdx='100-200ms_Oddball_25')
plotting.oddball10_si(which='MS', dest_dir_appdx='../plots/oddball_si', fname_appdx='100-200ms_Oddball_25vsMS')

# plotting.onset_offset_response(dest_dir_appdx='../plots/onset_offset_chnls')
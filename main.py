# coding: utf-8

import numpy as np
import warnings

from  MUA_utility import fetch, slice_data
import MUA_constants as const
from preprocessing import compress_CSVs
import plotting

def process_data():
    from MUA_cycle_dirs import MUA_analyzeAllFiles
    warnings.filterwarnings('ignore')
    MUA_analyzeAllFiles(delete_artifact_trials=True)
    warnings.filterwarnings('default')

# process_data()
# compress_CSVs() 

# plotting.firingrate_heatmaps('noisy', subtr_noise='paradigm_wise')
# plotting.firingrate_heatmaps('noise_subtr', 'dev_alone_C1C2')
# plotting.firingrate_noise_timeline('noisyy')
# plotting.firingrate_noise_timeline(fname_prefix='DA_subtr', subtr_noise='deviant_alone')
plotting.firingrate_noise_timeline(fname_prefix='parad_subtr', subtr_noise='paradigm_wise')

# plotting.plot_si('initial')
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
    MUA_analyzeAllFiles()
    warnings.filterwarnings('default')

# process_data()
# compress_CSVs() 
plotting.firingrate_heatmaps('noisy')
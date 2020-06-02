import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import glob

from MUA_utility import fetch, slice_data

# Investigating the differences between standard/deviant/pre_dev/post_dev
mouse = 'mGE85_31.07.2019_O10C1.mcd'

def preprocess_O10(p):
    data = fetch(p, ['mGE85'], ['O10C1'])
                        
    # for frates_fn in glob.glob(f'../output/output/*{mouse}/*_FiringRates.csv'):
    for key, frates in slice_data(data, firingrate=True).items():

        print(frates)
        plt.imshow(frates)
        plt.title(key)
        plt.show()




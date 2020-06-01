import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import glob

# Investigating the differences between standard/deviant/pre_dev/post_dev
mouse = 'mGE85_31.07.2019_O10C1.mcd'

def preprocess_O10(p):
    stim_types = ('Standard', 'Predeviant', 'Postdeviant', 'Deviant')
    data = dict().fromkeys(stim_types, None)
                        
    for frates_fn in glob.glob(f'../output/output/*{mouse}/*_FiringRates.csv'):
        which = [key for key in data.keys() if key in frates_fn][0]
        data[which] = pd.read_csv(frates_fn, index_col=0)
        # data[which] = data[which].reindex(p['id'][0]-1)
        print(data[which])
        
        plt.imshow(data[which])
        plt.title(which)
        plt.show()




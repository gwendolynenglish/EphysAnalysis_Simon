import pandas as pd
import numpy as np
import glob
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap

from MUA_utility import fetch, slice_data
from MUA_constants import ALL_MICE, ALL_PARADIGMS, ALL_STIMTYPES, PARAD_FULL

"""Investigate firingrates for different paradigm between 4 different mice"""
def explore_paradigms(p):
    def plot_paradigm(parad):
        data = fetch(p, paradigms=[parad])

        fig, axes = plt.subplots(4,4, sharex=True, sharey=True, figsize=(13,13))
        fig.subplots_adjust(hspace=.06, wspace=.03, right=.98, top=.86, left=.1, bottom=.07)
        # [ax.spines[where].set_visible(False) for ax in axes.flatten() for where in ax.spines]
        [ax.tick_params(bottom=False, left=False, labelbottom=False, labelleft=False) for ax in axes.flatten()]
        parad_full = PARAD_FULL[parad]
        fig.suptitle(f'{parad_full}- mean firing rates across 4 mice', size=14)
        plt.cm.get_cmap('gnuplot').set_gamma(.8)

        for mouse, i in zip(ALL_MICE, range(4)):
            mouse_dat = slice_data(data,[mouse], firingrate=True).items()
            axes[i,0].set_ylabel(mouse+'\nchannels', size=12, rotation=0, ha='right',
                    va='center')

            for (key, frates), j in zip(mouse_dat, range(4)):
                im = axes[i,j].imshow(frates, cmap='gnuplot', aspect='auto', extent=[-52.5, 202.5, -.5, 31.5],
                                    vmin=0, vmax=500)
                axes[i,j].vlines(0, -.5, 31.5, color='#ffffff', alpha=.6, linewidth=1)
                
                if i == 0:
                    axes[i,j].set_title(key[key.rfind('-')+1:])
                elif i == 3:
                    axes[i,j].tick_params(bottom=True, labelbottom=True)
                    axes[i,j].set_xlabel('ms')
                if (i == 0) and (j == 0):
                    axes[i,j].set_xlim((-52.5, 202.5))
                    axes[i,j].set_xticks([-50, 0, 80, 160])

                    # colorbar and legend
                    at = (0.77, .95, .2, .012,)
                    cb = fig.colorbar(im, cax=fig.add_axes(at), orientation='horizontal')
                    cb.set_label('Mean Firing Rate in 5ms frame', size=12)
                    cb.ax.get_xaxis().set_label_position('top')
        return fig

    for parad in ALL_PARADIGMS:
        fig = plot_paradigm(parad)
        path = p['outputPath']
        plt.savefig(f'{path}/../plots/allMice_{parad}_allStimTypes.png')
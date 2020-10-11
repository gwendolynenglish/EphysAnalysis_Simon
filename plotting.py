# Gwendolyn English 16.09.2019
# Functions for plotting LFPs  
################################################################################

################################################################################
# Inputs to the script are generally:
# 1. Evoked LFPs/MUAs 
# 2. Parameter dictionary 
################################################################################

# Import required Packages 
import numpy as np
import pandas as pd
import scipy.stats as ss
import pickle 
import matplotlib.mlab as mlab 
import matplotlib 
from matplotlib.pyplot import *
from statistics import * 

import os
import shutil
from glob import glob
from matplotlib import pyplot as plt
from PIL import Image

import MUA_constants as const
from shutil import move

# matplotlib.use('agg')
from analysisFunctionsLFP import down_sample_1D

################################################################################
def plot_evoked_shank(shankdata, outputpath, trigger, time, shank, shankmap): 
                      
    #Calculate boundaries
    max_evoked = np.max(shankdata)
    y_dist = 1.5 * max_evoked

    #Initialize figure
    fig = matplotlib.pyplot.figure()
    font = {'family': 'serif', 'color': 'black', 'weight': 'medium', 'size': 12}

    #Set y-ticks and labels 
    y_tickmarks = []
    y_labels = [] 
    for channel in np.arange(np.shape(shankdata)[0]):
        y_tickmarks.append(channel*y_dist)
        y_labels.append(str(np.shape(shankdata)[0] - channel))
        #y_labels.append(str(shankdata[channel]))

    yticks(y_tickmarks, shankmap[::-1] + 1)

    #Loop over channels in shank data
    for channel in np.arange(np.shape(shankdata)[0]):
        #plot(time, shankdata[channel] + (y_dist * channel))
        xdat = shankdata[int(np.shape(shankdata)[0] -channel-1)] + (y_dist*channel)
        plot(time, xdat)

    xlabel('Evoked Field Potentials', fontdict = font)
    ylabel('Channel', fontdict = font)

    savefig(outputpath + '/Evoked_shank_' + str(shank) + '_' + trigger + '.png', 
            format = 'png')
    savefig(outputpath + '/Evoked_shank_' + str(shank) + '_' + trigger + '.svg', 
            format = 'svg')
    close(fig)

################################################################################
def plot_evoked_channel(channelData, outputpathFolder, triggerFile, channelFile):    
    #Initialize figure
    fig = matplotlib.pyplot.figure()
    font = {'family': 'serif', 'color': 'black', 'weight': 'medium', 'size': 12}
    
    xlabel('Time (ms)', fontdict = font)
    ylabel('Evoked Response (uV)', fontdict = font) 
    xticks([0, len(channelData)-1],[const.P['evoked_pre']*-1000, const.P['evoked_post']*1000])
    title('Evoked Response ' + channelFile[:-4], fontdict = font) 
    
    channelData = channelData * 1000000
    plot(channelData)
    fname = outputpathFolder + '/Average_Evoked_Response_' + triggerFile[:-4] + \
            '_' + channelFile[:-4] + '.png'
    savefig(fname, format = 'png')
    close(fig)
    
################################################################################
def plot_raster(channelData, outputpath, artifact_trials):    
    
    #Initialize figure
    fig = matplotlib.pyplot.figure()
    font = {'family': 'serif', 'color': 'black', 'weight': 'medium', 'size': 12}
     
    if const.P['evoked_pre'] == 0.05 and const.P['evoked_post'] == 0.2: 
        xticks([-50,0,50,100,150,200])
        
    xlabel('Time (ms)', fontdict = font)
    ylabel('Stimulus', fontdict = font) 
    
    eventplot(channelData, colors = 'black', lineoffsets=1, linelengths=1 , 
              orientation = 'horizontal')
    
    if artifact_trials is not None:
        move(outputpath, f'{outputpath[:-4]}_PreArtifRemovel.png')
    savefig(outputpath, format = 'png') 
    close(fig)
    
################################################################################
def plot_PSTH(tsData, outputpath):    
    to = const.P['evoked_post'] * 1000 + const.P['psth_binsize']
    bins = np.arange(-const.P['evoked_pre'] *1000, to, const.P['psth_binsize'])
    histogram = np.histogram(tsData, bins) 
    
    #Initialize figure
    fig = matplotlib.pyplot.figure()
    font = {'family': 'serif', 'color': 'black', 'weight': 'medium', 'size': 12}
    
    xlabel('Time (ms)', fontdict = font)
    ylabel('Spike Count', fontdict = font) 
    
    hist(tsData[tsData !=0], bins)
    savefig(outputpath, format = 'png')
    close(fig)
    
################################################################################
def plot_waveforms(waveforms, outputpath):    

    #Initialize figure
    fig = matplotlib.pyplot.figure()
    font = {'family': 'serif', 'color': 'black', 'weight': 'medium', 'size': 12}

    xlabel('Time (ms)', fontdict = font)
    onems = const.P['sample_rate']/1000
    xticks([0, onems, onems*2, onems*3, onems*4] , [-2,-1,0,1,2])
    ylabel('Spike Amplitude (uV)', fontdict = font) 
    
    for spike in range(len(waveforms)):
        plot(waveforms[spike] *1000000)

    savefig(outputpath, format = 'png')
    close(fig)
    
################################################################################
def plot_firing_rate(firingrates, outputpath):    
    
    #Initialize figure
    fig = matplotlib.pyplot.figure()
    fig, ax = subplots()
    font = {'family': 'serif', 'color': 'black', 'weight': 'medium', 'size': 12}
    
    ax.set_title('Firing Rate', fontdict = font) 
    ax.set_xlabel('Time (ms)', fontdict = font)
    ax.set_ylabel('Firing Rate (Hz)', fontdict = font) 
    ax.set_xticks([-50,0,200])
    
    frm, to = -const.P['evoked_pre'] *1000, const.P['evoked_post'] * 1000
    bins = np.arange(frm, to, const.P['psth_binsize'])
    plot(bins, firingrates)
    axvline(x=0, color = 'r')

    savefig(outputpath, format = 'png')
    close(fig)
    
################################################################################
def plot_PSTH_bestfit_gaussian(tsData, outputpath):    
    frm, to = -const.P['evoked_pre'] *1000, const.P['evoked_post'] * 1000 + const.P['psth_binsize']
    bins = np.arange(frm, to, const.P['psth_binsize'])
    histogram = np.histogram(tsData, bins) 
    
    # Calculate mean and standard deviation of all post-stimulus spikes in 
    # histogram distribution 
    #Removes negative entries and turns into 1D array 
    data_clipped = tsData[tsData>0]    
    mu, sigma = ss.norm.fit(data_clipped)
    y = mlab.normpdf(bins, mu, sigma)
    
    #Initialize figure
    fig = matplotlib.pyplot.figure()
    font = {'family': 'serif', 'color': 'black', 'weight': 'medium', 'size': 12}
    
    #Formatting axes and title
    xlabel('Time (ms)', fontdict = font)
    ylabel('Spike Density', fontdict = font) 
    tit = (r'$\mathrm{Density\ of\ Evoked\ Spikes\ with\ Gaussian\ '
           'Fit:}\ \mu=%.2f,\ \sigma=%.2f$' %(mu, sigma))
    title(tit, fontdict = font)
    
    #Plot histogram and fit 
    hist(tsData[tsData !=0], bins, density = True)
    plot(bins, y, '--r')
    savefig(outputpath, format = 'png')
    close(fig)
    
################################################################################
def plot_PSTH_bestfit_gamma(tsData, outputpath):
    frm = -const.P['evoked_pre'] *1000, 
    to = const.P['evoked_post'] * 1000 + const.P['psth_binsize']
    bins = np.arange(frm, to, const.P['psth_binsize'])
    histogram = np.histogram(tsData, bins) 
    
    #Remove values below 0 for gamma fit 
    #Removes negative entries and turns into 1D array 
    data_clipped = tsData[tsData>0]    
    param = ss.gamma.fit(data_clipped, floc = 0)
    y = ss.gamma.pdf(bins, *param)
    k, alpha, theta = ss.gamma.fit(data_clipped, floc = 0)
    
    #Calculate moments of gamma fit
    mu = k*theta                   #mean
    sigma = np.sqrt(k * theta**2)  #std 
    gamma = 2/np.sqrt(k)           #skew 
    
    #Initialize figure
    fig = matplotlib.pyplot.figure()
    font = {'family': 'serif', 'color': 'black', 'weight': 'medium', 'size': 12}
    
    #Formatting axes and title
    xlabel('Time (ms)', fontdict = font)
    ylabel('Spike Density', fontdict = font) 
    title(r'$\mathrm{Density\ of\ Evoked\ Spikes\ with\ Gamma\ '
          'Fit:}\ \mu=%.2f,\ \sigma=%.2f,\ \gamma=%.2f$' %(mu, sigma, gamma))
    
    #Plot histogram and fit 
    hist(tsData[tsData !=0], bins, density = True)
    plot(y, '--r')
    savefig(outputpath, format = 'png')
    close(fig)

################################################################################
def plot_CSD(CSD, outputpath, trigger, time, shank, shankmap):

    #Calculate boundaries
    max_CSD = np.max(CSD)
    y_dist = 1.5 * max_CSD

    #Initialize figure
    fig = matplotlib.pyplot.figure()
    font = {'family': 'serif', 'color': 'black', 'weight': 'medium', 'size': 12}

    #Set y-ticks and labels 
    y_tickmarks = []
    y_labels = [] 
    for channel in np.arange(np.shape(CSD)[0]):
        y_tickmarks.append(channel*y_dist)
        y_labels.append(str(np.shape(CSD)[0] - channel))

    yticks(y_tickmarks, shankmap[::-1] + 1)

    #Loop over channels in CSD 
    for channel in np.arange(np.shape(CSD)[0]):
        #plot(time, CSD[channel] + (y_dist * channel))
        plot(time, CSD[int(np.shape(CSD)[0] - channel)-1] + (y_dist * channel))

    xlabel('Current Source Density', fontdict = font)
    ylabel('Channel', fontdict = font)

    savefig(outputpath + '/CSD_shank_' + str(shank) + '_' + trigger + '.png', 
            format = 'png')
    savefig(outputpath + '/CSD_shank_' + str(shank) + '_' + trigger + '.svg', 
            format = 'svg')
    close(fig)    

################################################################################
def plot_CSD_heatmap(CSD, outputpath, trigger, time, shank, shankmap):
    
    #Initialize figure
    fig = matplotlib.pyplot.figure()
    font = {'family': 'serif', 'color': 'black', 'weight': 'medium', 'size': 12}
    y = np.linspace(np.shape(CSD)[0],1,np.shape(CSD)[0])

    min = round(np.min(CSD))

    #Set y-ticks and labels 
    y_tickmarks = []
    y_labels = []
    for channel in np.arange(np.shape(CSD)[0]):
        y_tickmarks.append(channel)
    xticks([0, len(time)-1], [time[0] *1000, round(time[-1]*1000)])
    channels = shankmap + 1 
    #yticks(y_tickmarks, channels[::-1])
    yticks(y_tickmarks, channels[::1])

    #imshow(CSD[::-1], cmap = 'jet', aspect = 'auto', vmin = min, vmax = -min) 
    imshow(CSD[::1], cmap = 'jet', aspect = 'auto', vmin = min, vmax = -min) 
    colorbar(label = r'$\mu A / mm^{3}$')
    title('Current Source Density', fontdict = font)
    xlabel('Time(ms)', fontdict = font)
    ylabel('Channel', fontdict = font)
    savefig(outputpath + '/CSD_Heatmap_shank_' +str(shank) +'_'+ trigger + '.svg',  
            format = 'svg')
    savefig(outputpath + '/CSD_Heatmap_shank_' +str(shank) +'_'+ trigger + '.png',  
            format = 'png')
    close(fig)
    
################################################################################
def plot_wavelet_heatmap(avg_coef, freq, outputpath, triggerFile, channelFile, 
                         input, p=const.P):
    
    #Initialize figure
    fig = matplotlib.pyplot.figure()
    num = int(p['sample_rate']*(p['evoked_pre'] +p['evoked_post']))
    time = np.linspace(-p['evoked_pre'], p['evoked_post'], num)
    ds_time = down_sample_1D(time, p['sample_rate'])

    pcolor(ds_time, freq, avg_coef, cmap = 'seismic') 
    
    #ylim = ([1, p['low_pass_freq']])
    ylim = ([1,120])
    ylabel('Frequency (Hz)')
    xlabel('Time (ms)')
    colorbar()
    
    savefig(outputpath + '/Wavelet_Transform_' + triggerFile[:-4] + '_' + \
            channelFile[:-4] + input + '.svg', format = 'svg')
    savefig(outputpath + '/Wavelet_Transform_' + triggerFile[:-4] + '_' + \
            channelFile[:-4] + input + '.png', format = 'png')
    close(fig)













    
# --SIMON--
from matplotlib import pyplot as plt
from matplotlib.patches import Patch
from collections import OrderedDict

from  MUA_utility import fetch, slice_data, compute_si
import MUA_constants as const

def covariance_artifact_heatmap(cov, fault_trials, filename):
    plt.figure()
    plt.imshow(cov, aspect='auto', vmin=const.ARTIFACT_TRIAL_COV_HM_MIN, 
               vmax=const.ARTIFACT_TRIAL_COV_HM_MAX)
    plt.colorbar()

    fault_trials_idx = np.arange(0, len(fault_trials))[fault_trials]
    plt.xticks(fault_trials_idx)
    plt.yticks(fault_trials_idx)
    plt.gca().set_xticklabels([])
    plt.gca().set_yticklabels([])
    plt.xlabel(f'Trials ({cov.shape[0]})')
    plt.ylabel(f'Trials ({cov.shape[0]})')
    plt.title(f'Covariance of trials over electrical channels\n{sum(fault_trials)} deleted')

    plt.savefig(filename)

################################################################################
"""Investigate firingrates for different paradigm between different mice.
`subtr_noise` should either be False or `dev_alone_C1C2` or `deviant_alone` 
(ie the method)."""
def firingrate_heatmaps(dest_dir_appdx='../', subtr_noise=False, grouping='paradigm_wise'):
    def plot_paradigm(parad):
        if 'C1' in parad or 'C2' in parad:
            dev = parad[-2:]
            std = 'C2' if dev == 'C1' else 'C1'
        
        if grouping == 'paradigm_wise':
            data = fetch(paradigms=[parad])
        elif grouping == 'whisker_wise':
            if parad != 'MS':
                dev_data = fetch(paradigms=[parad], stim_types=['Deviant'])
                std_parad = parad.replace(dev, std) 
                std_data = fetch(paradigms=[std_parad], stim_types=['Standard', 'Predeviant', 'Postdeviant'])
                data = {**std_data, **dev_data}
            else:
                data = fetch(paradigms=[parad])
                
        elif grouping == 'whisker_wise_reduced':
            dev_parads = [this_parad for this_parad in const.ALL_PARADIGMS if dev in this_parad]
            std_parads = [this_parad.replace(dev, std) for this_parad in dev_parads]
            if dev == 'C1':
                order = ('DAC1-Deviant', 'O10C2-Predeviant', 'O10C1-Deviant', 'O25C2-Predeviant', 
                        'O25C1-Deviant', 'MS-C1', 'DOC2-Standard', 'DOC1-Deviant', )
            else:
                order = ('DAC2-Deviant', 'O10C1-Predeviant', 'O10C2-Deviant', 'O25C1-Predeviant', 
                        'O25C2-Deviant', 'MS-C2', 'DOC1-Standard', 'DOC2-Deviant')

            dev_data = fetch(paradigms=dev_parads, stim_types=['Deviant'])
            std_data = fetch(paradigms=std_parads, stim_types=['Predeviant'])
            std_data.update(fetch(paradigms=['DO'+std], stim_types=['Standard']))
            ms_data = fetch(paradigms=['MS'], stim_types=[dev])
            data = {**std_data, **dev_data, **ms_data}
            data = OrderedDict({key: data[key] for ord_key in order for key in data.keys() if ord_key in key})

        if grouping != 'whisker_wise_reduced':
            args = {'ncols': 4}
            width = 13
        else:
            args = {'ncols': 8 + 4, 'gridspec_kw': {'width_ratios': [.1, .015, .1, .1, .015, .1, .1, .015, .1, .015, .1, .1],}}
            width = 20
        fig, axes = plt.subplots(4, **args, sharex=True, sharey=True, figsize=(width,13))
        fig.subplots_adjust(hspace=.06, wspace=.03, right=.98, top=.86, left=.1, bottom=.07)
        
        [ax.tick_params(bottom=False, left=False, labelbottom=False, labelleft=False) for ax in axes.flatten()]
        if grouping != 'whisker_wise_reduced':
            title = const.PARAD_FULL[parad] + '- mean firing rates across 4 mice'
        else:
            title = parad + '- mean firing rates across 4 mice'

        if subtr_noise:
            title += '\n\nNOISE SUBTRACTED'
        fig.suptitle(title, size=14)
        plt.cm.get_cmap('gnuplot').set_gamma(.8)

 
        for mouse, i in zip(const.ALL_MICE, range(4)):
            mouse_dat = slice_data(data, [mouse], firingrate=True, 
                                   frate_noise_subtraction=subtr_noise)
            axes[i,0].set_ylabel(mouse+'\nchannels', size=12, rotation=0, ha='right',
                    va='center')

            which_ax = 0
            for (key, frates), j in zip(mouse_dat.items(), range(args['ncols'])):
                
                im = axes[i,which_ax].imshow(frates, cmap='gnuplot', aspect='auto', extent=[-52.5, 202.5, -.5, 31.5],
                                    vmin=0, vmax=500)
                axes[i,which_ax].vlines(0, -.5, 31.5, color='#ffffff', alpha=.6, linewidth=1)
                
                if i == 0:
                    stim_t = key[key.rfind('-')+1:]
                    if grouping == 'paradigm_wise' or parad == 'MS':
                        print(stim_t)
                        axes[i,which_ax].set_title(stim_t)
                    elif grouping == 'whisker_wise'  :
                        axes[i,which_ax].set_title('C1 '+ stim_t)
                    elif grouping == 'whisker_wise_reduced':
                        pard_full = const.PARAD_FULL[key[key.find('-')+1:key.rfind('-')]][:-3]
                        title = f'{dev} {stim_t}\n{pard_full}'
                        if 'MS' not in key:
                            axes[i,which_ax].set_title(title)
                        else:
                            axes[i,which_ax].set_title(stim_t+'\nMany Standards')

                elif i == 3:
                    axes[i,which_ax].tick_params(bottom=True, labelbottom=True)
                    axes[i,which_ax].set_xlabel('ms')
                if (i == 0) and (which_ax == 0):
                    axes[i,which_ax].set_xlim((-52.5, 202.5))
                    axes[i,which_ax].set_xticks([-50, 0, 80, 160])

                    # colorbar and legend
                    at = (0.77, .95, .2, .012,)
                    cb = fig.colorbar(im, cax=fig.add_axes(at), orientation='horizontal')
                    cb.set_label('Mean Firing Rate in 5ms frame', size=12)
                    cb.ax.get_xaxis().set_label_position('top')
                
                which_ax += 1
                if grouping == 'whisker_wise_reduced' and which_ax in [1,4,7,9]:
                    axes[i, which_ax].set_visible(False)
                    which_ax += 1
        return fig
    
    paradigms = const.ALL_PARADIGMS
    if grouping == 'whisker_wise_reduced':
        paradigms = 'C1', 'C2'
    for parad in paradigms:
        fig = plot_paradigm(parad)
        plt.savefig(f'{const.P["outputPath"]}/{dest_dir_appdx}/{parad}_allStimTypes.png')
        plt.close(fig)

def firingrate_noise_timeline(fname_prefix='', subtr_noise=False):
    """Check background firingrates activity for experimental time line between 
    4 different mice and averaged stimulus types. `subtr_noise` should either 
    be False or `deviant_alone`, `paradigm_wise` (ie the method)."""

    data = fetch()

    ratio = {'width_ratios': [.1] *11 + [.28],}
            #  'height_ratios': [.15, .03, .8]}
    fig, axes = plt.subplots(4,12, sharex=False, sharey=False, figsize=(15,13), gridspec_kw=ratio)
    fig.subplots_adjust(hspace=.2, wspace=.03, right=.95, top=.86, left=.1, bottom=.07)

    title = 'noise over experimental paradigm timeline\n(mean over stimulus types)'
    if subtr_noise:
        title += '\n\nNoise subtracted: ' + subtr_noise
    fig.suptitle(title, size=14)
    plt.cm.get_cmap('gnuplot').set_gamma(.8)
    [ax.tick_params(bottom=False, left=False, labelbottom=False, labelleft=False) for ax in axes.flatten()]

    for mouse, i in zip(const.ALL_MICE, range(4)):
        mouse_dat = slice_data(data, [mouse], firingrate=True, 
                                frate_noise_subtraction=subtr_noise)
        axes[i,0].set_ylabel(mouse+'\nchannels', size=12, rotation=0, ha='right',
                            va='center')

        neg_frates = []
        for parad, j in zip(const.PARAD_ORDER[mouse], range(11)):
            parad_frates = [df for key, df in mouse_dat.items() if parad in key]
            
            neg_frates_counts = [(fr<0).sum().sum() for fr in parad_frates]
            neg_frates.append(sum(neg_frates_counts) /len(neg_frates_counts))
            
            frates = (sum(parad_frates)/ len(parad_frates)).astype(int)
            im = axes[i,j].imshow(frates, cmap='gnuplot', aspect='auto', 
                                  extent=[-52.5, 202.5, -.5, 31.5], vmin=0, 
                                  vmax=500)

            axes[i,j].set_title(parad, size=7, pad=2)
            if 'DA' in parad and subtr_noise == 'deviant_alone':
                axes[i,j].set_title('**'+parad+'**', size=9, pad=2)
                [axes[i,j].spines[where].set_color('yellow') for where in axes[i,j].spines]
                [axes[i,j].spines[where].set_linewidth(1.5) for where in axes[i,j].spines]

            if (i == 0) and (j == 0):
                axes[i,j].set_xlim((-52.5, 202.5))
                axes[i,j].set_xticks([-50, 0, 80, 160])

                # colorbar and legend
                at = (0.77, .95, .2, .012,)
                cb = fig.colorbar(im, cax=fig.add_axes(at), orientation='horizontal')
                cb.set_label('Mean Firing Rate in 5ms frame', size=12)
                cb.ax.get_xaxis().set_label_position('top')
        
        rel_neg_frates = np.array(neg_frates) /1600
        axes[i,11].bar(range(11), rel_neg_frates)
        axes[i,11].bar([12], rel_neg_frates.mean())

        # x axis        
        axes[i,11].tick_params(labelbottom=True, rotation=35, labelsize=6.5)
        axes[i,11].set_xticks(list(range(11))+[12])
        axes[i,11].set_xticklabels(const.PARAD_ORDER[mouse]+('Average',), clip_on=False, ha='right', y=.04)
        
        # y axis
        axes[i,11].tick_params(right=True, labelright=True)
        axes[i,11].set_ylim((0,1))
        yticks = np.arange(0,1.1,.1)
        axes[i,11].set_yticks(yticks)
        axes[i,11].set_yticklabels([str(yt) if yt in (0,1) else '' for yt in yticks])
        axes[i,11].yaxis.set_label_position("right")
        axes[i,11].yaxis.grid(True, alpha=.5)
        axes[i,11].set_axisbelow(True)

        if i == 1:
            axes[i,11].set_ylabel('Proportion of negative firing rates (of 32 channels x 50 time bins)', size=10)
    
    path = const.P['outputPath']
    plt.savefig(f'{path}/../plots/firingrates_lowthr/{fname_prefix}_firingrate_noise_over_time.png')



















def oddball10_si(dest_dir_appdx, fname_appdx, which='O10'):
    data = fetch(mouseids = ['mGE82', 'mGE83', 'mGE84', 'mGE85'], 
                 paradigms = [which+'C1', which+'C2', 'MS'] if not which == 'MS' else ['O25C1', 'O25C2', 'MS'], 
                 stim_types = ['Deviant', 'Predeviant', 'C1', 'C2'], collapse_ctx_chnls=True, 
                 collapse_th_chnls=True, drop_not_assigned_chnls=True)
    
    SIs, frates = compute_si(data, MS=which=='MS',  start=120, stop=160)
    # plt.hist(frates, bins=200)
    # print(frates)
    # print(SIs.isna().sum())
    # plt.show()
    SIs_mean = SIs.mean()

    fig, ax = plt.subplots(figsize=(8, 6))
    fig.subplots_adjust(top=.75, right=.82, left=.2, bottom=.15)

    [sp.set_visible(False) for sp in ax.spines.values()]
    ax.spines['left'].set_visible(True)
    ax.patch.set_facecolor('grey')
    ax.patch.set_alpha(.16)
    ax.hlines((0),0,23, color='black', linewidth=.5)
    ax.set_title(fname_appdx + ' SSA', pad=65)

    xt = [1,2, 6,7, 11,12, 16,17, 21,22]
    xt_mid = [1.5, 6.5, 11.5, 16.5, 21.5]
    ax.set_xlim((0,23))
    ax.set_xticks(xt)
    ax.set_xticklabels(['C2' if i%2 else 'C1' for i in range(len(xt))])

    ax.yaxis.grid(True, which='major')
    ax.set_ylim((-1.05,1.05))
    ax.set_ylabel('SSA index (SI)')
    ax.set_yticks(np.arange(-1, 1.001, .25))
    
    colors = [const.COLORS[col] for col in ['red', 'deep_blue', 'magenta','teal']]
    for (m_id, mouse_si), col in zip(SIs.iterrows(), colors):
        ax.scatter(xt, mouse_si, color=col, s=6, alpha=.5, label=m_id)

    regions = [const.REGIONS_EXT[reg] for reg in SIs_mean.index.unique(0)]
    [ax.annotate(reg, (x_m, 1.05), rotation=30) for reg, x_m in zip(regions, xt_mid)]
    ax.scatter(xt, SIs_mean, color='k', s=20, marker='x', label='Average')
    ax.legend(bbox_to_anchor=(1.001, 1.001), loc='upper left')

    lbl = f'avg. firingr. in 5ms < {const.SI_MIN_FRATE_5MS}\n# excluded mice:'
    ax.annotate(lbl, (-7.3,-1.5), annotation_clip=False, fontsize=9)
    [ax.annotate(n, (x_t-.3, -1.5), annotation_clip=False, fontsize=9) for n, x_t in zip(SIs.isna().sum().values, xt)]
    SIs.isna().sum()

    f = f'{const.P["outputPath"]}/{dest_dir_appdx}/SSA_indices_{fname_appdx}.png'
    fig.savefig(f)


def ssa_correlation(dest_dir_appdx, fname_appdx, which='O10', post_stim=False):
    data = fetch(mouseids = ['mGE82', 'mGE83', 'mGE84', 'mGE85'], 
                 paradigms = [which+'C1', which+'C2', 'MS'] if not which == 'MS' else ['O25C1', 'O25C2', 'MS'], 
                 stim_types = ['Deviant', 'Predeviant', 'C1', 'C2'], collapse_ctx_chnls=True, 
                 collapse_th_chnls=True, drop_not_assigned_chnls=True)
    
    SIs, frates = compute_si(data, MS=which=='MS', start=5, stop=20)
    
    if post_stim:
        SIs_post, frates_post = compute_si(data, MS=which=='MS', start=100, stop=200)
        SIs_post.columns = pd.MultiIndex.from_tuples([(region_whisker[0]+'_lateSI', region_whisker[1]) for region_whisker in SIs_post.columns])
        frates_post.columns = pd.MultiIndex.from_tuples([(region_whisker[0]+'_lateSI', region_whisker[1]) for region_whisker in frates_post.columns])
        order = SIs.columns.unique(0).tolist() + SIs_post.columns.unique(0).tolist()

        SIs = pd.concat([SIs, SIs_post], axis=1)
        frates = pd.concat([frates, frates_post], axis=1)

        SIs = SIs.stack(level=1).reindex(order, axis=1)
        frates = frates.stack(level=1).reindex(order, axis=1)
    else:
        SIs = SIs.stack(level=1).reindex(['Th', 'G', 'SG', 'IG', 'dIG'], axis=1)
        frates = frates.stack(level=1).reindex(['Th', 'G', 'SG', 'IG', 'dIG'], axis=1)

    corr = np.corrcoef(SIs.T)

    fig, ax = plt.subplots(figsize=(6, 6))
    im = ax.imshow(corr, aspect='auto', vmin=-1, vmax=1, cmap='RdBu_r')
    
    # colorbar and legend
    at = (0.77, .95, .2, .012,)
    cb = fig.colorbar(im, cax=fig.add_axes(at), orientation='horizontal')
    cb.set_label('Spearman\'s r', size=12)
    cb.ax.get_xaxis().set_label_position('top')

    ax.set_title(f'SSA index correlation {fname_appdx}')
    ax.set_xticks(np.arange(SIs.shape[1]))
    ax.set_xticklabels(SIs.columns, fontsize=10, rotation=45)
    ax.set_yticks(np.arange(SIs.shape[1]))
    ax.set_yticklabels(SIs.columns, fontsize=10, rotation=45, rotation_mode='anchor')
    f = f'{const.P["outputPath"]}/{dest_dir_appdx}/SSA_corr_heatmap_{fname_appdx}.png'
    fig.savefig(f)

    for comp_reg, comp_dat in SIs.iteritems():
        for i, (reg, region_dat) in enumerate(SIs.iteritems()):
            if reg == comp_reg:
                continue
            fig, ax = plt.subplots(figsize=(6, 6))

            [sp.set_visible(False) for sp in ax.spines.values()]
            ax.patch.set_facecolor('grey')
            ax.patch.set_alpha(.16)
            ax.hlines((0),-1,1, color='black', linewidth=.5)
            ax.vlines((0),-1,1, color='black', linewidth=.5)

            ax.set_xlim(-.75,1.05)
            ax.set_ylim(-.3,1.05)

            ax.set_xlabel('SSA index '+const.REGIONS_EXT[comp_reg])
            ax.set_ylabel('SSA index '+const.REGIONS_EXT[reg])
            
            ax.scatter(comp_dat, region_dat,s=5, color='k')
            # if comp_reg == 'Th':
            #     [ax.annotate(frates.loc[idx, comp_reg], (comp_dat[idx], region_dat[idx]), size=7) for idx in frates[reg].index]
            # [ax.annotate(frates.loc[idx, reg], (comp_dat[idx], region_dat[idx]), size=7, ha='right', va='top') for idx in frates[comp_reg].index]
            [ax.annotate('-'.join(idx), (comp_dat[idx], region_dat[idx]), size=7, ha='right', 
                         va='bottom' if 'C1' in idx else 'top') for idx in frates[comp_reg].index]

            r = ss.linregress(comp_dat, region_dat)
            ax.plot((-1,0,1), (r.intercept-r.slope, r.intercept, r.slope+r.intercept), 
                   color=const.REGION_CMAP[reg], label=f'{reg} p-value: {r.pvalue:.2f}')
            plt.legend(loc='lower left')

            for idx in region_dat.index:
                reg_d = region_dat.drop(idx)
                comp_d = comp_dat.drop(idx)
                
                plt.scatter(comp_d, reg_d, color=const.REGION_CMAP[reg],s=3)
                r = ss.linregress(comp_d, reg_d)

                ax.plot((-1,0,1), (r.intercept-r.slope, r.intercept, r.slope+r.intercept), 
                        linestyle=(0, (5, 10)), linewidth=.8, color=const.REGION_CMAP[reg], label=f'{reg} p-value: {r.pvalue:.2f}')

            f = f'{const.P["outputPath"]}/{dest_dir_appdx}/SSA_corr_{comp_reg}-{reg}_{fname_appdx}.png'
            fig.savefig(f)


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
def firingrate_heatmaps(dest_dir_appdx='../', subtr_noise=False, all_stimuli=False,
                        grouping='whisker_wise', chnls_to_regions=False):
    """Generate a panel of firing rate heatmaps for all mice. The specific split
    regarding the paradigm is controlled by `grouping`. By default, this equates
    to `whisker_wise` which will make a plot for each paradigm but importantly,
    pooled such that only one whisker is stimulated. This means for example for 
    the panel of O10C1, the deviant of O10C1 is plotted plus the Standard in 
    O10C2. This way, the effect of the paradigm on one whsker, here C1, is 
    ideally visualized. When grouping equals `whisker_wise_reduced`, the 
    Predeviant and Postdeviant are cut, but one big panel showing all paradigms
    for each whisker is produced. This is gives a high level overview of 
    differences between paradigms. Lastly, grouping may be passed as 
    `paradigm_wise` which simply takes an experiment and plots the
    Deviant and Standard found in there. For O10C1, the Deviant column would 
    corresbond to C1 stimulation, the Standard columns to C2 stimulation. This 
    is an old option that is not really that useful.  `subtr_noise` may be 
    passed as False (default), 'paradigm_wise' (recommended) or 'deviant_alone'.
    When `chnls_to_regions` is passed as True (defaults to False), the channel
    mapping file is used to collapse channels to regions. The heatmaps then have
    5 rows instead of 32. `all_stimuli` will include postdevient and standard 
    stimuli, by default, only predevient and devient are plotted (2 columns).
    Currently only implemented for grouping=`whisker_wise`.
    Plots will be saved as usual at P["outputPath"]/dest_dir_appdx/* """
    def plot_paradigm(parad):
        print(parad)
        if 'C1' in parad or 'C2' in parad:
            dev = parad[-2:]
            std = 'C2' if dev == 'C1' else 'C1'
        
        to_regions = ['collapse_ctx_chnls', 'collapse_th_chnls', 'drop_not_assigned_chnls']
        to_regions = dict.fromkeys(to_regions, True if chnls_to_regions else False)
        
        if grouping == 'paradigm_wise':
            data = fetch(paradigms=[parad])
        elif grouping == 'whisker_wise':
            if parad != 'MS':
                dev_data = fetch(paradigms=[parad], stim_types=['Deviant'], **to_regions)
                
                std_parad = parad.replace(dev, std) 
                stim_types = ['Predeviant', 'UniquePredeviant']
                if all_stimuli:
                    stim_types.extend(['Standard', 'Postdeviant', 'UniquePostdeviant'])
                std_data = fetch(paradigms=[std_parad], stim_types=stim_types, **to_regions)
                if std_parad in ['DOC1', 'DOC2']:
                    std_data = fetch(paradigms=[std_parad], stim_types=['Standard'], **to_regions)

                data = {**std_data, **dev_data}
            else:
                data = fetch(paradigms=[parad], **to_regions)
                
        elif grouping == 'whisker_wise_reduced':
            dev_parads = [this_parad for this_parad in const.ALL_PARADIGMS if dev in this_parad]
            std_parads = [this_parad.replace(dev, std) for this_parad in dev_parads]
            if dev == 'C1':
                order = ('DAC1-Deviant', 'O10C2-Predeviant', 'O10C1-Deviant', 'O25C2-UniquePredeviant', 
                        'O25C1-Deviant', 'MS-C1', 'DOC2-Standard', 'DOC1-Deviant', )
            else:
                order = ('DAC2-Deviant', 'O10C1-Predeviant', 'O10C2-Deviant', 'O25C1-UniquePredeviant', 
                        'O25C2-Deviant', 'MS-C2', 'DOC1-Standard', 'DOC2-Deviant')

            dev_data = fetch(paradigms=dev_parads, stim_types=['Deviant'])
            std_data = fetch(paradigms=std_parads, stim_types=['Predeviant', 'UniquePredeviant'])
            std_data.update(fetch(paradigms=['DO'+std], stim_types=['Standard']))
            ms_data = fetch(paradigms=['MS'], stim_types=[dev])
            data = {**std_data, **dev_data, **ms_data}
            data = OrderedDict({key: data[key] for ord_key in order for key in data.keys() if ord_key in key})

        if grouping != 'whisker_wise_reduced':
            if parad != 'MS' and not all_stimuli:
                args = {'ncols': 2}
                width = 7
            else:
                args = {'ncols': 4}
                width = 13

        else:
            args = {'ncols': 8 + 4, 'gridspec_kw': {'width_ratios': [.1, .015, .1, .1, .015, .1, .1, .015, .1, .015, .1, .1],}}
            width = 20
        fig, axes = plt.subplots(4, **args, sharex=True, sharey=True, figsize=(width,13))
        fig.subplots_adjust(hspace=.06, wspace=.03, right=.98, top=.86, left=.13, bottom=.07)
        
        [ax.tick_params(bottom=False, left=False, labelbottom=False, labelleft=False) for ax in axes.flatten()]
        if grouping != 'whisker_wise_reduced':
            title = const.PARAD_FULL[parad] + '- mean firing rates across 4 mice'
        else:
            title = parad + '- mean firing rates across 4 mice'

        if subtr_noise:
            title += '\nNOISE SUBTRACTED'
        fig.suptitle(title, size=14)
        plt.cm.get_cmap('gnuplot').set_gamma(.8)

        print()
        print()
        print()
        for mouse, i in zip(const.ALL_MICE, range(4)):
            mouse_dat = slice_data(data, [mouse], firingrate=True, 
                                   frate_noise_subtraction=subtr_noise)
            if chnls_to_regions:
                axes[i,0].tick_params(left=True, labelleft=True)
                lbls = list(mouse_dat.values())[0].index[::-1]
                axes[i,0].set_yticks(np.arange(6.4/2, 32, 6.4))
                axes[i,0].set_yticklabels(lbls)

            axes[i,0].set_ylabel(mouse+'\nchannels', size=12, rotation=0, ha='right',
                    va='center', x=-50)

            which_ax = 0
            for (key, frates), j in zip(mouse_dat.items(), range(args['ncols'])):
                
                im = axes[i,which_ax].imshow(frates, cmap='gnuplot', aspect='auto', extent=[-52.5, 202.5, -.5, 31.5],
                                    vmin=0, vmax=500)
                axes[i,which_ax].vlines(0, -.5, 31.5, color='w', alpha=.6, linewidth=1)
                if i == 0:
                    stim_t = key[key.rfind('-')+1:]
                    col = 'k' if stim_t != 'Deviant' else const.COLORS['deep_red']
                    if grouping == 'paradigm_wise' or parad == 'MS':
                        axes[i,which_ax].set_title(stim_t, color=col)
                    elif grouping == 'whisker_wise'  :
                        axes[i,which_ax].set_title(f'{parad[-2:]} {stim_t}', color=col)
                    elif grouping == 'whisker_wise_reduced':
                        pard_full = const.PARAD_FULL[key[key.find('-')+1:key.rfind('-')]][:-3]
                        title = f'{dev} {stim_t}\n{pard_full}'
                        if 'MS' not in key:
                            axes[i,which_ax].set_title(title, color=col)
                        else:
                            axes[i,which_ax].set_title(stim_t+'\nMany Standards', color=col)

                elif i == 3:
                    axes[i,which_ax].tick_params(bottom=True, labelbottom=True)
                    axes[i,which_ax].set_xlabel('ms')
                if (i == 0) and (which_ax == 0):
                    axes[i,which_ax].set_xlim((-52.5, 202.5))
                    axes[i,which_ax].set_xticks([-50, 0, 80, 160])

                    # colorbar and legend
                    at = (.58, .9, 2.5/width, .012,)
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


def firingrate_noise_timeline(dest_dir_appdx='../', fname_postfix='', subtr_noise=False):
    """Check background firingrates activity for experimental time line between 
    4 different mice and averaged stimulus types. `subtr_noise` should either 
    be False or `deviant_alone`, `paradigm_wise` (ie the method)."""

    data = fetch()

    ratio = {'width_ratios': [.1] *11 + [.28],}
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
    
    plt.savefig(f'{const.P["outputPath"]}/{dest_dir_appdx}/firingrate_over_time{fname_postfix}.{const.PLOT_FORMAT}')

def oddball_si(dest_dir_appdx, which='O10', compare_with_MS=False, start=5, stop=20):
    """ Plot the SSA index of a specific paradigm. The paradigm is passed to 
    `which`, it may be `O10`, `O25` or `O25U`. By default the SSA index is 
    calculated by comparing the devient response to the predevient (standard). 
    The response is defined as the average firing rate in a given time interval
    post stimulus. This interval is passed by `start` (default 5) and `stop` 
    (default 20). The SSA index may also be calculated using the many standards
    paradigm. Instead of comparing with the standard presentation of the whisker,
    the stimulation of the whisker within the MS paradigm is used. This is done
    by passing `compare_with_MS` as True. Details on the compuation of SSA 
    indecis can be found in the doctring of the MUA_utility.py function 
    `compute_si()`. Besides the main plot, a histogram of respones is saved to
    check the general magnitute of responses. All plots are saved at 
    P["outputPath"]/dest_dir_appdx/*.
    """
    versusMS = 'versusMS' if compare_with_MS else ''
    data = fetch(mouseids = ['mGE82', 'mGE83', 'mGE84', 'mGE85'], 
                 paradigms = [which+'C1', which+'C2', 'MS'], 
                 stim_types = ['Deviant', 'Predeviant', 'UniquePredeviant', 'C1', 'C2'], 
                 collapse_ctx_chnls=True, collapse_th_chnls=True, 
                 drop_not_assigned_chnls=True)

    SIs, _, SI_raw_values = compute_si(data, MS=compare_with_MS,  start=start, stop=stop)
    SIs_mean = SIs.mean()
    
    fig, (ax_left, ax_right) = plt.subplots(ncols=2, figsize=(7, 5))
    [ax_left.hist(SI_raw_values[m_id].values.flatten(), label=m_id, alpha=.8, bins=20, range=(0,10), color=const.GENERAL_CMAP[m_id]) for m_id in SI_raw_values.columns.unique(0)]
    [ax_right.hist(SI_raw_values.values.flatten(), bins=40)]
    ax_left.set_ylim(ax_left.get_ylim())
    ax_left.set_ylabel('counts')
    ax_left.set_xlabel(f'avg. firingrate {start}-{stop} ms')
    ax_right.set_xlabel(f'avg. firingrate {start}-{stop} ms')
    ax_left.set_title('low responses')
    ax_right.set_title('all responses')
    ax_left.legend()
    ax_left.vlines(const.SI_MIN_FRATE_5MS, -5, 50, linestyle='dashed')
    ax_left.annotate('cut off', (const.SI_MIN_FRATE_5MS+.1, ax_left.get_ylim()[1]*.9))
    fig.suptitle(f'SSA {which} {versusMS} {start}-{stop} ms')
    
    versusMS = 'versusMS' if compare_with_MS else ''
    f = f'{const.P["outputPath"]}/{dest_dir_appdx}/responses_hist_{which}_{versusMS}_{start}-{stop}ms.{const.PLOT_FORMAT}'
    fig.savefig(f)

    fig, ax = plt.subplots(figsize=(8, 6))
    fig.subplots_adjust(top=.75, right=.82, left=.2, bottom=.15)

    [sp.set_visible(False) for sp in ax.spines.values()]
    ax.spines['left'].set_visible(True)
    ax.patch.set_facecolor('grey')
    ax.patch.set_alpha(.16)
    ax.hlines((0),0,23, color='black', linewidth=.5)
    ax.set_title(f'SSA {which} {versusMS} {start}-{stop} ms', pad=65)

    xt = [1,2, 6,7, 11,12, 16,17, 21,22]
    xt_mid = [1.5, 6.5, 11.5, 16.5, 21.5]
    ax.set_xlim((0,23))
    ax.set_xticks(xt)
    ax.set_xticklabels(['C2' if i%2 else 'C1' for i in range(len(xt))])

    ax.yaxis.grid(True, which='major')
    ax.set_ylim((-1.05,1.05))
    ax.set_ylabel('SSA index (SI)')
    ax.set_yticks(np.arange(-1, 1.001, .25))
    
    for m_id, mouse_si in SIs.iterrows():
        ax.scatter(xt, mouse_si, color=const.GENERAL_CMAP[m_id], s=8, alpha=.9, label=m_id)

    regions = [const.REGIONS_EXT[reg] for reg in SIs_mean.index.unique(0)]
    [ax.annotate(reg, (x_m, 1.05), rotation=30) for reg, x_m in zip(regions, xt_mid)]
    ax.scatter(xt, SIs_mean, color='k', s=20, alpha=.7, marker='D', label='Average')
    ax.legend(bbox_to_anchor=(1.001, 1.001), loc='upper left')

    lbl = f'avg. firingr. in 5ms < {const.SI_MIN_FRATE_5MS}\n# excluded mice:'
    ax.annotate(lbl, (-7.3,-1.5), annotation_clip=False, fontsize=9)
    [ax.annotate(n, (x_t-.3, -1.5), annotation_clip=False, fontsize=9) for n, x_t in zip(SIs.isna().sum().values, xt)]
    SIs.isna().sum()

    versusMS = 'versusMS' if compare_with_MS else ''
    f = f'{const.P["outputPath"]}/{dest_dir_appdx}/SSA_index_{which}_{versusMS}_{start}-{stop}ms.{const.PLOT_FORMAT}'
    fig.savefig(f)

def ssa_correlation(dest_dir_appdx, which='O10', start=5, stop=20, post_stim=False):
    """Plot the correlation of SSA indices. Each mouse has two SIs, one for C1,
    one for C2. With 4 mice, we get a vector of 8 datapoints that are then 
    correlated in different regions. If one of those 8 is NA, it is excluded 
    from the correlation. At least 4 values must be present for the correlation 
    of two regions. A heatmap is drawn showing pearson's rho coefficient. To
    check if the correlation looks reasonable, for each combination of regions,
    a plot is drawn that shows the SIs and a fitted regression line. The 
    paradigm is passed in `which`. `start` and `stop` also work as mentioned 
    before. `post_stim` may simply be passed as True to include late responses 
    that are hard fixed in the code to late_start=100 and late_stop=200. This 
    essentially extends the regions domain from G, IG .. to G, IG, ... late_G, 
    late IG.... Often for late responses the cut off response constant 
    should be lowered. Therefore, post_stim can also be a number that will 
    be subtracted from the value defined in MUA_constants.py SI_MIN_FRATE_5MS. 
    Plots are saved at P["MUA_ouput"]/dest_dir_appdx/*
    """
    data = fetch(mouseids = ['mGE82', 'mGE83', 'mGE84', 'mGE85'], 
                 paradigms = [which+'C1', which+'C2', 'MS'] if not which == 'MS' else ['O25C1', 'O25C2', 'MS'], 
                 stim_types = ['Deviant', 'Predeviant', 'UniquePredeviant', 'C1', 'C2'], collapse_ctx_chnls=True, 
                 collapse_th_chnls=True, drop_not_assigned_chnls=True)
    
    SIs, frates, _ = compute_si(data, MS=which=='MS', start=start, stop=stop)
    
    if post_stim:
        late_start=100
        late_stop=200
        print(const.SI_MIN_FRATE_5MS)
        if post_stim is not True:
            const.SI_MIN_FRATE_5MS -= post_stim
        print(const.SI_MIN_FRATE_5MS)
        SIs_post, frates_post, _ = compute_si(data, MS=which=='MS', start=late_start, stop=late_stop)
        SIs_post.columns = pd.MultiIndex.from_tuples([(region_whisker[0]+'_lateSI', region_whisker[1]) for region_whisker in SIs_post.columns])
        frates_post.columns = pd.MultiIndex.from_tuples([(region_whisker[0]+'_lateSI', region_whisker[1]) for region_whisker in frates_post.columns])
        order = SIs.columns.unique(0).tolist() + SIs_post.columns.unique(0).tolist()

        SIs = pd.concat([SIs, SIs_post], axis=1)
        frates = pd.concat([frates, frates_post], axis=1)

        SIs = SIs.stack(level=1).reindex(order, axis=1)
        frates = frates.stack(level=1).reindex(order, axis=1)
    else:
        SIs = SIs.stack(level=1).reindex(['VPM', 'G', 'SG', 'IG', 'dIG'], axis=1)
        frates = frates.stack(level=1).reindex(['VPM', 'G', 'SG', 'IG', 'dIG'], axis=1)

    p_values = {}
    late = 'late' if post_stim else ''
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
            [ax.annotate('-'.join(idx), (comp_dat[idx], region_dat[idx]), size=7, ha='right', 
                         va='bottom' if 'C1' in idx else 'top') for idx in frates[comp_reg].index]

            notna = comp_dat.notna().values & region_dat.notna()
            if notna.sum() <=4:
                p_values[f'{comp_reg}-{reg}'] = 'NaN'
            else:
                r = ss.linregress(comp_dat[notna], region_dat[notna])
                ax.plot((-1,0,1), (r.intercept-r.slope, r.intercept, r.slope+r.intercept), 
                        color=const.REGION_CMAP[reg], label=f'{reg} p-value: {r.pvalue:.2f}')
                p_values[f'{comp_reg}-{reg}'] = r.pvalue
                plt.legend(loc='lower left')

            f = f'{const.P["outputPath"]}/{dest_dir_appdx}/SSA_corr_{comp_reg}-{reg}_{which}_{start}_{stop}ms_{late}.{const.PLOT_FORMAT}'
            fig.savefig(f)
    
    print(SIs.to_string())
    SIs.T[SIs.notna().sum()<=4] = np.nan
    print(SIs.to_string())
    corr = SIs.corr()

    figsize = (8,8) if not post_stim else (11,11)
    fig, ax = plt.subplots(figsize=figsize)
    fig.subplots_adjust(left=.1, bottom=.1, top=.75, right=.75)
    im = ax.imshow(corr, aspect='auto', vmin=-1, vmax=1, cmap='RdBu_r')
    
    for row, reg in enumerate(SIs.columns):
        for col, reg_nd in enumerate(SIs.columns):
            if reg == reg_nd:
                continue
            pval = p_values[f'{reg}-{reg_nd}']
            pval = f'p={pval:.3f}' if type(pval) is not str else 'NaN'
            ax.annotate(pval, (row-.35, col), fontsize=8)
    
    # colorbar and legend
    at = (0.77, .95, .2, .012,)
    cb = fig.colorbar(im, cax=fig.add_axes(at), orientation='horizontal')
    cb.set_label('Pearson\'s r', size=12)
    cb.ax.get_xaxis().set_label_position('top')

    ax.set_title(f'SSA index correlation {which} {start}-{stop}ms')
    ax.set_xticks(np.arange(SIs.shape[1]))
    ax.set_xticklabels(SIs.columns, fontsize=10, rotation=45)
    ax.set_yticks(np.arange(SIs.shape[1]))
    ax.set_yticklabels(SIs.columns, fontsize=10, rotation=45, rotation_mode='anchor')

    n_smples = [f'{reg} n={nsmples}' for reg, nsmples in  SIs.notna().sum().iteritems()]
    ax.annotate('\n'.join(n_smples), (.77, .6), annotation_clip=False, xycoords='figure fraction')

    f = f'{const.P["outputPath"]}/{dest_dir_appdx}/SSA_corr_heatmap_{which}_{start}_{stop}ms_{late}.{const.PLOT_FORMAT}'
    fig.savefig(f)
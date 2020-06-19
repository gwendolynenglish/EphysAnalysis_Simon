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
def plot_raster(channelData, outputpath):    
    
    #Initialize figure
    fig = matplotlib.pyplot.figure()
    font = {'family': 'serif', 'color': 'black', 'weight': 'medium', 'size': 12}
     
    if const.P['evoked_pre'] == 0.05 and const.P['evoked_post'] == 0.2: 
        xticks([-50,0,50,100,150,200])
        
    xlabel('Time (ms)', fontdict = font)
    ylabel('Stimulus', fontdict = font) 
    
    eventplot(channelData, colors = 'black', lineoffsets=1, linelengths=1 , 
              orientation = 'horizontal')
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

from  MUA_utility import fetch, slice_data, compute_si
import MUA_constants as const

################################################################################
"""Investigate firingrates for different paradigm between 4 different mice.
`subtr_noise` should either be False or `dev_alone_C1C2` (ie the method)."""
def firingrate_heatmaps(fname_prefix='', subtr_noise=False):
    def plot_paradigm(parad):
        data = fetch(paradigms=[parad])

        fig, axes = plt.subplots(4,4, sharex=True, sharey=True, figsize=(13,13))
        fig.subplots_adjust(hspace=.06, wspace=.03, right=.98, top=.86, left=.1, bottom=.07)
        
        [ax.tick_params(bottom=False, left=False, labelbottom=False, labelleft=False) for ax in axes.flatten()]
        title = const.PARAD_FULL[parad] + '- mean firing rates across 4 mice'
        if subtr_noise:
            title += '\n\nNOISE SUBTRACTED'
        fig.suptitle(title, size=14)
        plt.cm.get_cmap('gnuplot').set_gamma(.8)

        for mouse, i in zip(const.ALL_MICE, range(4)):
            mouse_dat = slice_data(data,[mouse], firingrate=True, 
                                   frate_noise_subtraction=subtr_noise)
            axes[i,0].set_ylabel(mouse+'\nchannels', size=12, rotation=0, ha='right',
                    va='center')

            for (key, frates), j in zip(mouse_dat.items(), range(4)):
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

    for parad in const.ALL_PARADIGMS:
        fig = plot_paradigm(parad)
        path = const.P['outputPath']
        plt.savefig(f'{path}/../plots/firingrates/{fname_prefix}_allMice_{parad}_allStimTypes.png')
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
    plt.savefig(f'{path}/../plots/firingrates/{fname_prefix}_firingrate_noise_over_time.png')


def plot_si(fname_prefix):
    data = fetch(paradigms=['O10C1', 'O10C2'])
    
    SI_values = compute_si(data)

    std = [val for key, val in SI_values.items() if 'Standard' in key]
    pred = [val for key, val in SI_values.items() if 'Predeviant' in key]
    postd = [val for key, val in SI_values.items() if 'Postdeviant' in key]


    fig, ax = plt.subplots()
    ax.set_xticks([1,2,3])
    ax.set_xticklabels(['Standard', 'Predeviant', 'Postdeviant'])
    ax.set_ylabel('SI index')
    ax.set_title('Oddball 10% - C1,C2 mean')
    ax.set_ylim((0,1))

    for i in range(4):
        ax.scatter([1,2,3], [std[i], pred[i], postd[i]])

    plt.show()



def make_CSD_summary_plots(lfp_output_appdx, dest_dir_appdx):
    """ p['outputPath']/`dest_dir_appdx --> output_folder
        p['outputPath']/`LFP_output_appdx --> LFP_analysis output folder
    This function copies and creates a collection of CSD relevent plots and
    makes one panel to compare them easily. For each mouse|paradigm, a subfolder
    is created in `dest_dir` :=  p['outputPath']/`dest_dir_appdx and populated
     with the CSD heatmap & CSD lineplot (copied over), a scatter plot 
    showing the amplitude and timestamp of the first LFP post stim, and lastly, 
    an average firingrate heatmap for each channel 100ms post stim. Finally, A 
    table is printed out that contains the assigned cortical region. This table 
    is initialized empty at `dest_dir`/ctx_mapping.csv. The idea is to go over 
    these summary panels one by one and fill the channel-to-cortical region 
    table along the way. ctx_mapping.csv should be filled wth values SG, G, IG, 
    dIG. When this function is rerun with the table populated, the table on the 
    summary panels updates to show the assigned cortical region and therefore 
    serves as a handy documentation how the channels where assigned. Fiannly,
    the summary plots are copied into the all_panels directory. 
    """

    # cp the CSD heatmap and the line plot over to the mouse|paradigm 
    # target directory where all the plots are collected
    def copy_over_CSD_plots(CSD_heatmap, CSD_lineplot):
        plot_images = []
        for plot_file in (CSD_heatmap, CSD_lineplot):
            plot_images.append(f'{dest_dir}/{mouse_parad_dir}/{os.path.basename(plot_file)}')
            shutil.copyfile(plot_file, plot_images[-1])
        # return the save locations
        return plot_images
    
    # plot the timestamp vs the amplitude of the first negative LFP for each channel
    def draw_neg_peak_plot(CSD_lfp_avg_csv):
        plt.figure(figsize=(8,8))
        plt.subplots_adjust(bottom=.1, top=.9, left=0.1, right=.95)
        plt.xlim((0, 200))
        plt.xlabel('time [ms]', fontsize=14)
        plt.ylabel(' first max neg ampl. [uV]', fontsize=14)
        plt.title(f'Layer 4 identification\n{parad_order}: {mouse_parad_dir}', fontsize=14)
        
        lfp_avg = pd.read_csv(CSD_lfp_avg_csv, index_col=0)
        time_stamps = lfp_avg.Peak_neg_ts
        ampl = lfp_avg.Peak_neg_uV
        plt.scatter(time_stamps, ampl, s=5)
        [plt.annotate(int(lfp_avg.Electrode[i]), (time_stamps[i], ampl[i]), fontsize=11) 
         for i in range(lfp_avg.shape[0])]
        
        plot_filename = (f'{dest_dir}/{mouse_parad_dir}/negative_peak_poststim_plot.png')
        plt.savefig(plot_filename, tight_layout=True)
        plt.close()
        return plot_filename

    # draw a heatmap visualizing the binned firingrate as in the main function above 
    def draw_avg_frate_plot(firingrates_csv):
        plt.figure(figsize=(4,8))
        plt.subplots_adjust(bottom=0.01, top=.9, left=0.1, right=.98)
        plt.tick_params(bottom=False, labelbottom=False, top=True, labeltop=True) 
        
        plt.title('firing rates post stim')
        plt.xlabel('time post stim [ms]', labelpad=-10)
        plt.xticks([0, 100])
        plt.gca().xaxis.set_label_position('top')
        plt.yticks(np.arange(32), [f'{lbl:.0f}' for lbl in const.P['id'][0]][::-1])
        
        frates = pd.read_csv(firingrates_csv, index_col=0)
        frates = frates.iloc(1)[10:30]
        plt.imshow(frates, cmap='gnuplot', aspect='auto', extent=[-2.5, 102.5, -.5, 31.5],
                                    )
        plt.hlines(np.arange(-.5, 31.5), -2.5, 102.5, color='w', alpha=.3)
        
        plot_file = f'{dest_dir}/{mouse_parad_dir}/frate_barh.png'
        plt.savefig(plot_file, tight_layout=True)
        plt.close()
        return plot_file

    # the table is simply a plt.barh plot that is colored according to a color
    # mapping defined below 
    def draw_ctx_mapping(ctx_mapping_csv):
        if not os.path.exists(ctx_mapping_csv):
            cols = [f'{mid}-{parad}' for mid in const.PARAD_ORDER.keys() 
                    for parad in const.PARAD_ORDER[mid]]
            ctx_df = pd.DataFrame('not_assigned', columns=cols, 
                                  index=const.P['id'][0])
            ctx_df.to_csv(ctx_mapping_csv)
        ctx_map_df = pd.read_csv(ctx_mapping_csv, index_col=0)
    
        plt.figure(figsize=(2,8))
        plt.subplots_adjust(bottom=.01, top=.9, left=0.02, right=.83)
        plt.tick_params(bottom=False, left=False, labelleft=False, right=True, labelright=True)
        plt.xlim((0,1))
        plt.ylim((31.5,-.5))
        plt.yticks(np.arange(32), [f'{lbl:.0f}' for lbl in const.P['id'][0]])
        
        mapping = ctx_map_df[f'{mouse}-{paradigm}']
        # define how each layer is colored 
        ctx_color_map = {'not_assigned': 'w', 'SG': '#42d4f4', 'G': '#e6194B', 
                        'IG': '#bfef45', 'dIG': '#aaffc3'}
        cols = [ctx_color_map[region] for region in mapping]
        plt.barh(np.arange(32), 1, height=1, edgecolor='k', color=cols, alpha=.6)
        [plt.annotate(mapping.values[i], (.2, i+.2), fontsize=10) for i in range(len(mapping))]
        
        plot_file = f'{dest_dir}/{mouse_parad_dir}/ctx_map_table.png'
        plt.savefig(plot_file, tight_layout=True)
        plt.close()
        return plot_file

    # the composition function uses the PIL library to make one big white panel,
    # then pastes the collected images at specific coordinates 
    def make_image_comp(plot_images):
        imgs = [Image.open(img) for img in plot_images]
        height = imgs[0].height + imgs[2].height
        width = imgs[0].width + imgs[1].width + imgs[4].width
        final_img = Image.new('RGB', (width, height))
        
        final_img.paste(imgs[0], (0, 0))    # upper left
        final_img.paste(imgs[1], (imgs[0].width, 0))    # upper right
        final_img.paste(imgs[2], (0, imgs[0].height))   # lower left
        final_img.paste(imgs[3], (imgs[2].width, imgs[0].height))   # lower middle
        final_img.paste(imgs[4], (imgs[2].width+imgs[3].width, imgs[0].height)) # lower right
        return final_img

    
    dest_dir = const.P['outputPath'] + dest_dir_appdx
    lfp_output = const.P['outputPath'] + lfp_output_appdx
    
    all_panels = []
    # get the last part of the LFP_output directories, ie the mouse-date-paradigm.mcd
    mouse_paradigms = [os.path.basename(path) for path in glob(f'{lfp_output}/mGE*')]
    for mouse_parad_dir in mouse_paradigms:
        mouse = mouse_parad_dir[:mouse_parad_dir.find('_')]
        paradigm = mouse_parad_dir[mouse_parad_dir.rfind('_')+1:-4]
        parad_order = [i+1 for i, parad in enumerate(const.PARAD_ORDER[mouse]) 
                       if parad == paradigm][0]
        
        # create a mouse|paradigm specific directory
        os.makedirs(f'{dest_dir}/{mouse_parad_dir}', exist_ok=True)

        # define all the files to get (input)
        CSD_heatmap = f'{lfp_output}/{mouse_parad_dir}/CSD_Heatmap_shank_0_Triggers_Deviant.dat.png'
        CSD_lineplot = f'{lfp_output}/{mouse_parad_dir}/Evoked_shank_0_Triggers_Deviant.dat.png'
        CSD_lfp_avg_csv = f'{lfp_output}/{mouse_parad_dir}/Triggers_Deviant_LFPAverages.csv'
        firingrates_csv = f'{lfp_output}/../MUA_output/{mouse_parad_dir}/Triggers_Deviant_FiringRates.csv'

        # all metrics always depend on the deviant stimulus, if MS get C1,
        # if DOC get Standard 
        if 'MS' in mouse_parad_dir or 'DOC' in mouse_parad_dir:
            instead = 'C1' if 'MS' in mouse_parad_dir else 'Standard'
            CSD_heatmap = CSD_heatmap.replace('Deviant', instead)
            CSD_lineplot = CSD_lineplot.replace('Deviant', instead)
            CSD_lfp_avg_csv = CSD_lfp_avg_csv.replace('Deviant', instead)
            firingrates_csv = firingrates_csv.replace('Deviant', instead)

        # filepaths for ctx_mapping.csv and the final panel
        ctx_mapping_csv = f'{dest_dir}/ctx_mapping.csv'
        # add paradaigm order to the final filename to see possible effects over time
        CSD_composition = f'{dest_dir}/{mouse_parad_dir}/CSD_summay_plot_{mouse}_{parad_order}_{paradigm}.png'

        # call the functions defined on top one by one, eg copy CSD_heatmap and
        # CSD_lineplot, draw 2 new plots and the table, always get the filenames 
        # back and save in `plot_images`
        plot_images = copy_over_CSD_plots(CSD_heatmap, CSD_lineplot)
        plot_images.append(draw_neg_peak_plot(CSD_lfp_avg_csv))
        plot_images.append(draw_avg_frate_plot(firingrates_csv))
        plot_images.append(draw_ctx_mapping(ctx_mapping_csv))
        
        # make the final panel
        final_img = make_image_comp(plot_images)
        final_img.save(CSD_composition)
        all_panels.append(CSD_composition)
    
    # make a dir with all the summry plots symlinked in
    os.makedirs(f'{dest_dir}/all_panels', exist_ok=True)
    [shutil.copyfile(src, f'{dest_dir}/all_panels/{os.path.basename(src)}') 
     for src in all_panels]

def plot_evoked_lfp(lfp_output_appdx, dest_dir_appdx):
    path = const.P['outputPath'] + lfp_output_appdx + '/mGE82_24.07.2019_DAC1.mcd'
    mouse = 'mGE82'
    paradigm = 'DAC1'
    ctx_mapping_csv = const.P['outputPath'] + '/../CSD/ctx_mapping.csv'
    data = pd.read_csv(path+'/Triggers_Deviant_LFPAverages.csv', index_col=0).iloc(1)[7:]
    # firingrates_csv = f'{lfp_output}/../MUA_output/{mouse_parad_dir}/Triggers_Deviant_FiringRates.csv'

    thal_data = data.iloc[-10:, :] *1_000_000
    ctx_data = data.iloc[:10, :] *1_000_000
    x_time = data.columns.values.astype(float)
    # print(df)
    # y = df.loc[0]
    # print(y)
    # plt.plot(y)
    # plt.show()

    # init figure
    fig = plt.figure(figsize=(15,12))
    gs = fig.add_gridspec(7, 3, width_ratios=[.5, .35, .15], 
                          height_ratios=[.25, .025, .15, .15, .15, .15, .15], 
                          hspace=0, wspace=.12, right=.93, top=.95, left=.07,
                          bottom=.05)
    lfp_ax_ctx = fig.add_subplot(gs[0, 0])
    frate_ax = fig.add_subplot(gs[0, 1])
    assig_ax = fig.add_subplot(gs[0, 2])
    lfp_axl = [fig.add_subplot(gs[i, 0]) for i in range(2, 7)]
    lfp_axr = [fig.add_subplot(gs[i, 1:]) for i in range(2, 7)]
    all_axs = lfp_axl + lfp_axr + [lfp_ax_ctx, frate_ax, assig_ax]
    lfp_axs = lfp_axl + lfp_axr + [lfp_ax_ctx]

    [ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False) 
     for ax in all_axs]

    for ax in lfp_axs:
        [sp.set_visible(False) for sp in ax.spines.values()]
        ax.patch.set_facecolor('grey')
        ax.patch.set_alpha(.16)
        ax.xaxis.grid(True, which='major')
        ax.yaxis.grid(True, which='major')
        ax.hlines(0,-5,20, color='grey')
        
        xleft, xright = -0.05, 0.25
        ax.set_xlim((xleft, xright))
        xticks = (xleft, 0, 0.05, 0.1, 0.15, 0.2, xright)
        ax.set_xticks(xticks)

        ybot, ytop = -25, 25
        ax.set_ylim((ybot, ytop))
        ax.set_yticks((-15, 0, 15))
    
    [(ax.tick_params(labelbottom=True), ax.set_xticklabels(xticks, size=8), ax.set_xlabel('[s]', size=8)) 
     for ax in (lfp_axl[-1], lfp_axr[-1])]
    [(ax.tick_params(labelright=True), ax.set_yticklabels((-15, 0, 15), size=8))  for ax in lfp_axr]
    lfp_axr[2].annotate('[uV]', (.24, 5), clip_on=False)
    

    for (chnl, dat), ax in zip(thal_data.iterrows(), lfp_axl+lfp_axr):
        ax.set_ylabel(chnl+1, fontsize=14, rotation=0, ha='right', va='center', labelpad=1)
        ax.plot(x_time, dat.values, clip_on=False, color='#469990')

    
    lfp_ax_ctx.set_ylim(-175, 100)
    yticks = np.arange(-175, 150, 15)
    lfp_ax_ctx.set_yticks(yticks)
    [(ax.tick_params(labelright=True), )  for ax in [lfp_ax_ctx]]

    
    SG = ctx_data.iloc[1:4, :].mean(0)
    G = ctx_data.iloc[4:6, :].mean(0)
    IG = ctx_data.iloc[6:12, :].mean(0)

    lfp_ax_ctx.plot(x_time, SG, label='SG', clip_on=False)
    lfp_ax_ctx.plot(x_time, G, label='G', clip_on=False)
    lfp_ax_ctx.plot(x_time, IG, label='IG', clip_on=False)
    lfp_ax_ctx.legend()

    ctx_map_df = pd.read_csv(ctx_mapping_csv, index_col=0)
    
    # plt.figure(figsize=(2,8))
    # plt.subplots_adjust(bottom=.01, top=.9, left=0.02, right=.83)
    assig_ax.tick_params(bottom=False, left=False, labelleft=False, right=True, labelright=True)
    assig_ax.set_xlim((0,1))
    assig_ax.set_ylim((31.5, 21.5))
    assig_ax.set_yticks(np.arange(31.5, 21.5, -1), [f'{lbl:.0f}' for lbl in const.P['id'][0]])
    
    mapping = ctx_map_df[f'{mouse}-{paradigm}']
    # define how each layer is colored 
    ctx_color_map = {'not_assigned': 'w', 'SG': '#42d4f4', 'G': '#e6194B', 
                    'IG': '#bfef45', 'dIG': '#aaffc3'}
    cols = [ctx_color_map[region] for region in mapping]
    assig_ax.barh(np.arange(32), 1, height=1, edgecolor='k', color=cols, alpha=.6)
    [assig_ax.annotate(mapping.values[i], (.2, i+.2), fontsize=10) for i in range(len(mapping))]
    
    # plot_file = f'{dest_dir}/{mouse_parad_dir}/ctx_map_table.png'
    # plt.savefig(plot_file, tight_layout=True)
    # plt.close()
    # return plot_file



        # def draw_avg_frate_plot(firingrates_csv):
    
    
    # plt.figure(figsize=(4,8))
    # plt.subplots_adjust(bottom=0.01, top=.9, left=0.1, right=.98)
    frate_ax.tick_params(bottom=False, labelbottom=False, top=True, labeltop=True) 
    
    # plt.title('firing rates post stim')
    # plt.xlabel('time post stim [ms]', labelpad=-10)
    # frate_ax.set_xticks([0, 100])
    frate_ax.xaxis.set_label_position('top')
    frate_ax.set_yticks(np.arange(32), [f'{lbl:.0f}' for lbl in const.P['id'][0]][::-1])
    
    frates = pd.read_csv(firingrates_csv, index_col=0)
    frates = frates.iloc(1)[10:30]
    plt.imshow(frates, cmap='gnuplot', aspect='auto', extent=[-2.5, 102.5, -.5, 31.5],
                                )
    # plt.hlines(np.arange(-.5, 31.5), -2.5, 102.5, color='w', alpha=.3)
        


        # plot_file = f'{dest_dir}/{mouse_parad_dir}/frate_barh.png'
        # plt.savefig(plot_file, tight_layout=True)
        # plt.close()
        # return plot_file

    plt.show()
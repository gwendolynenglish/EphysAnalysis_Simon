
##########################################################################################################
# Gwendolyn English 16.09.2019
# Functions for plotting LFPs  
##########################################################################################################

##########################################################################################################
# Inputs to the script are generally:
# 1. Evoked LFPs/MUAs 
# 2. Parameter dictionary 
##########################################################################################################

##########################################################################################################
# Import required Packages 
import numpy as np
import scipy.stats as ss
import pickle 
import matplotlib.mlab as mlab 
import matplotlib 
matplotlib.use('agg')
from matplotlib.pyplot import *
from statistics import * 
from analysisFunctionsLFP import * 

##########################################################################################################
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
        plot(time, shankdata[int(np.shape(shankdata)[0] - channel - 1)] + (y_dist * channel))

    xlabel('Evoked Field Potentials', fontdict = font)
    ylabel('Channel', fontdict = font)

    savefig(outputpath + '/Evoked_shank_' + str(shank) + '_' + trigger + '.png', format = 'png')
    savefig(outputpath + '/Evoked_shank_' + str(shank) + '_' + trigger + '.svg', format = 'svg')
    close(fig)

##########################################################################################################
def plot_evoked_channel(p, channelData, outputpathFolder, triggerFile, channelFile):    
    
    #Initialize figure
    fig = matplotlib.pyplot.figure()
    font = {'family': 'serif', 'color': 'black', 'weight': 'medium', 'size': 12}
    
    xlabel('Time (ms)', fontdict = font)
    ylabel('Evoked Response (uV)', fontdict = font) 
    xticks([0,  len(channelData)-1],[p['evoked_pre']*-1000, p['evoked_post']*1000])
    title('Evoked Response ' + channelFile[:-4], fontdict = font) 
    
    channelData = channelData * 1000000
    plot(channelData)
    savefig(outputpathFolder + '/Average_Evoked_Response_' + triggerFile[:-4] + '_' + channelFile[:-4] + '.png', format = 'png')
    close(fig)
    
##########################################################################################################
def plot_raster(p, channelData, outputpath):    
    
    #Initialize figure
    fig = matplotlib.pyplot.figure()
    font = {'family': 'serif', 'color': 'black', 'weight': 'medium', 'size': 12}
     
    if p['evoked_pre'] == 0.05 and p['evoked_post'] == 0.2: 
        xticks([-50,0,50,100,150,200])
        
    xlabel('Time (ms)', fontdict = font)
    ylabel('Stimulus', fontdict = font) 
    
    eventplot(channelData, colors = 'black', lineoffsets=1, linelengths=1 , orientation = 'horizontal')
    savefig(outputpath, format = 'png') 
    close(fig)
    
##########################################################################################################
def plot_PSTH(p, tsData, outputpath):    
    bins = np.arange(-p['evoked_pre'] *1000, p['evoked_post'] * 1000 + p['psth_binsize'], p['psth_binsize'])
    histogram = np.histogram(tsData, bins) 
    
    #Initialize figure
    fig = matplotlib.pyplot.figure()
    font = {'family': 'serif', 'color': 'black', 'weight': 'medium', 'size': 12}
    
    xlabel('Time (ms)', fontdict = font)
    ylabel('Spike Count', fontdict = font) 
    
    hist(tsData[tsData !=0], bins)
    savefig(outputpath, format = 'png')
    close(fig)
    
##########################################################################################################
def plot_waveforms(p, waveforms, outputpath):    

    #Initialize figure
    fig = matplotlib.pyplot.figure()
    font = {'family': 'serif', 'color': 'black', 'weight': 'medium', 'size': 12}

    xlabel('Time (ms)', fontdict = font)
    onems = p['sample_rate']/1000
    xticks([0, onems, onems*2, onems*3, onems*4] , [-2,-1,0,1,2])
    ylabel('Spike Amplitude (uV)', fontdict = font) 
    
    for spike in range(len(waveforms)):
        plot(waveforms[spike] *1000000)

    savefig(outputpath, format = 'png')
    close(fig)
    
##########################################################################################################
def plot_firing_rate(p, firingrates, outputpath):    
    
    #Initialize figure
    fig = matplotlib.pyplot.figure()
    font = {'family': 'serif', 'color': 'black', 'weight': 'medium', 'size': 12}
    
    xlabel('Time (ms)', fontdict = font)
    ylabel('Firing Rate (Hz)', fontdict = font) 
    #xticks([0,  len(firingrates)-1],[p['evoked_pre']*-1000, p['evoked_post']*1000])
    xticks([0, np.int((p['evoked_pre'] / p['evoked_post']) * len(firingrates)), len(firingrates)-1], \
           [p['evoked_pre']*-1000, 0, p['evoked_post']*1000])
    title('Firing Rate', fontdict = font) 
    
    #Identify t=0 and plot 
    stim_onset = np.int((p['evoked_pre'] / p['evoked_post']) * len(firingrates))
    axvline(x=stim_onset, color = 'r')
    
    plot(firingrates)
    savefig(outputpath, format = 'png')
    close(fig)
    
##########################################################################################################
def plot_PSTH_bestfit_gaussian(p, tsData, outputpath):    
    bins = np.arange(-p['evoked_pre'] *1000, p['evoked_post'] * 1000 + p['psth_binsize'], p['psth_binsize'])
    histogram = np.histogram(tsData, bins) 
    
    #Calculate mean and standard deviation of all post-stimulus spikes in histogram distribution 
    data_clipped = tsData[tsData>0]    #Removes negative entries and turns into 1D array 
    mu, sigma = ss.norm.fit(data_clipped)
    y = mlab.normpdf(bins, mu, sigma)
    
    #Initialize figure
    fig = matplotlib.pyplot.figure()
    font = {'family': 'serif', 'color': 'black', 'weight': 'medium', 'size': 12}
    
    #Formatting axes and title
    xlabel('Time (ms)', fontdict = font)
    ylabel('Spike Density', fontdict = font) 
    title(r'$\mathrm{Density\ of\ Evoked\ Spikes\ with\ Gaussian\ Fit:}\ \mu=%.2f,\ \sigma=%.2f$' %(mu, sigma), fontdict = font)
    
    #Plot histogram and fit 
    hist(tsData[tsData !=0], bins, density = True)
    plot(bins, y, '--r')
    savefig(outputpath, format = 'png')
    close(fig)
    
##########################################################################################################
def plot_PSTH_bestfit_gamma(p, tsData, outputpath):    
    bins = np.arange(-p['evoked_pre'] *1000, p['evoked_post'] * 1000 + p['psth_binsize'], p['psth_binsize'])
    histogram = np.histogram(tsData, bins) 
    
    #Remove values below 0 for gamma fit 
    data_clipped = tsData[tsData>0]    #Removes negative entries and turns into 1D array 
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
    title(r'$\mathrm{Density\ of\ Evoked\ Spikes\ with\ Gamma\ Fit:}\ \mu=%.2f,\ \sigma=%.2f,\ \gamma=%.2f$' %(mu, sigma, gamma))
    
    #Plot histogram and fit 
    hist(tsData[tsData !=0], bins, density = True)
    plot(y, '--r')
    savefig(outputpath, format = 'png')
    close(fig)

##########################################################################################################    
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

    savefig(outputpath + '/CSD_shank_' + str(shank) + '_' + trigger + '.png', format = 'png')
    savefig(outputpath + '/CSD_shank_' + str(shank) + '_' + trigger + '.svg', format = 'svg')
    close(fig)    

##########################################################################################################       
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
    savefig(outputpath + '/CSD_Heatmap_shank_' + str(shank) + '_' + trigger + '.svg',  format = 'svg')
    savefig(outputpath + '/CSD_Heatmap_shank_' + str(shank) + '_' + trigger + '.png',  format = 'png')
    close(fig)
    
##########################################################################################################       
def plot_wavelet_heatmap(avg_coef, freq, p, outputpath, triggerFile, channelFile, input):
    
    #Initialize figure
    fig = matplotlib.pyplot.figure()
    time = np.linspace(-p['evoked_pre'], p['evoked_post'], p['sample_rate']*(p['evoked_pre'] +p['evoked_post']))
    ds_time = down_sample_1D(time, p['sample_rate'])

    pcolor(ds_time, freq, avg_coef, cmap = 'seismic') 
    
    #ylim = ([1, p['low_pass_freq']])
    ylim = ([1,120])
    ylabel('Frequency (Hz)')
    xlabel('Time (ms)')
    colorbar()
    
    savefig(outputpath + '/Wavelet_Transform_' + triggerFile[:-4] + '_' + channelFile[:-4] + input + '.svg', format = 'svg')
    savefig(outputpath + '/Wavelet_Transform_' + triggerFile[:-4] + '_' + channelFile[:-4] + input + '.png', format = 'png')
    close(fig)
          
    
    
# coding: utf-8

import numpy as np
import warnings
import concurrent.futures
import os


from  MUA_utility import fetch, slice_data
import MUA_constants as const
from preprocessing import compress_CSVs
import plotting


"""Process data """
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



"""General explorative / summarizing  plots"""
# plotting.firingrate_heatmaps(dest_dir_appdx='../plots/frates_whiskerwise', subtr_noise='paradigm_wise', grouping='whisker_wise')
# plotting.firingrate_heatmaps(dest_dir_appdx='../plots/frates_whiskerwise_reduced', subtr_noise='paradigm_wise', grouping='whisker_wise_reduced')
# plotting.firingrate_heatmaps('noisy', False)
# plotting.firingrate_noise_timeline('noisyy')
# plotting.firingrate_noise_timeline(fname_prefix='parad_subtr', subtr_noise='paradigm_wise')



"""Mapping cortical and thalamic channels"""
# plotting.make_CSD_summary_plots(lfp_output_appdx='/../LFP_output', dest_dir_appdx='/../CSD_lowthr')
# plotting.plot_time_to_first(dest_dir=f'{const.PROJ_DIR}/output/time_to_first')
# plotting.plot_evoked_lfp(dest_dir_appdx='../thalamus_mapping_lowthr', 
#                          anatomy_dir=f'{const.PROJ_DIR}/metadata/Histology_pngs',
                        #  ts_dir=f'{const.PROJ_DIR}/output/time_to_first')



"""SSA and SSA correlation"""
# start, stop = 5, 20
# plotting.oddball10_si(which='O10', dest_dir_appdx='../plots/oddball_si', fname_appdx=f'{start}-{stop}ms_Oddball_10')
# plotting.oddball10_si(which='O25', dest_dir_appdx='../plots/oddball_si', fname_appdx=f'{start}-{stop}ms_Oddball_25')
# plotting.oddball10_si(which='MS', dest_dir_appdx='../plots/oddball_si', fname_appdx=f'{start}-{stop}ms_Oddball_25vsMS')
# plotting.ssa_correlation(which='O10', dest_dir_appdx='../plots/si_corr_5-20ms_010_late', fname_appdx='', post_stim=True)
# plotting.ssa_correlation(which='O10', dest_dir_appdx='../plots/si_corr_5-20ms_010_check', fname_appdx='', post_stim=False)



"""Onset offset"""
plotting.onset_offset_response(dest_dir_appdx='../plots/onset_offset_chnls')
plotting.onset_offset_response(dest_dir_appdx='../plots/onset_offset_regions', single_channels=False)
# plotting.find_onset_offsets('../plots/on_off_rasters', './onset_offset_scores_diag.csv')

"""
2 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O10C1.mcd/Raster_NegativeSpikes_Triggers_Standard_ElectrodeChannel_07.png

2 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O10C1.mcd/Raster_NegativeSpikes_Triggers_Standard_ElectrodeChannel_08.png

2 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O10C1.mcd/Raster_NegativeSpikes_Triggers_Standard_ElectrodeChannel_09.png

2 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O10C1.mcd/Raster_NegativeSpikes_Triggers_Predeviant_ElectrodeChannel_06.png

2 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O10C1.mcd/Raster_NegativeSpikes_Triggers_Predeviant_ElectrodeChannel_07.png

2 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O10C1.mcd/Raster_NegativeSpikes_Triggers_Predeviant_ElectrodeChannel_08.png

2 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O10C1.mcd/Raster_NegativeSpikes_Triggers_Predeviant_ElectrodeChannel_09.png

2 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O10C1.mcd/Raster_NegativeSpikes_Triggers_Postdeviant_ElectrodeChannel_07.png

2 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O10C1.mcd/Raster_NegativeSpikes_Triggers_Postdeviant_ElectrodeChannel_08.png

2 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O10C2.mcd/Raster_NegativeSpikes_Triggers_Standard_ElectrodeChannel_31.png

2 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE85_31.07.2019_O10C2.mcd/Raster_NegativeSpikes_Triggers_Standard_ElectrodeChannel_31.png

2 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE85_31.07.2019_O10C2.mcd/Raster_NegativeSpikes_Triggers_Standard_ElectrodeChannel_32.png

1 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE85_31.07.2019_O10C2.mcd/Raster_NegativeSpikes_Triggers_Predeviant_ElectrodeChannel_32.png

2 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE85_31.07.2019_O10C2.mcd/Raster_NegativeSpikes_Triggers_Postdeviant_ElectrodeChannel_31.png

1 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE85_31.07.2019_O10C2.mcd/Raster_NegativeSpikes_Triggers_Postdeviant_ElectrodeChannel_32.png

1 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O25C1.mcd/Raster_NegativeSpikes_Triggers_Standard_ElectrodeChannel_06.png

1 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O25C1.mcd/Raster_NegativeSpikes_Triggers_Standard_ElectrodeChannel_07.png

1 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O25C1.mcd/Raster_NegativeSpikes_Triggers_Standard_ElectrodeChannel_08.png

2 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O25C1.mcd/Raster_NegativeSpikes_Triggers_Standard_ElectrodeChannel_09.png

1 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O25C1.mcd/Raster_NegativeSpikes_Triggers_UniquePredeviant_ElectrodeChannel_06.png

1 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O25C1.mcd/Raster_NegativeSpikes_Triggers_UniquePredeviant_ElectrodeChannel_07.png

1 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O25C1.mcd/Raster_NegativeSpikes_Triggers_UniquePredeviant_ElectrodeChannel_08.png

1 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O25C1.mcd/Raster_NegativeSpikes_Triggers_UniquePredeviant_ElectrodeChannel_09.png

2 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O25C1.mcd/Raster_NegativeSpikes_Triggers_UniquePostdeviant_ElectrodeChannel_06.png

1 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O25C1.mcd/Raster_NegativeSpikes_Triggers_UniquePostdeviant_ElectrodeChannel_07.png

1 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O25C1.mcd/Raster_NegativeSpikes_Triggers_UniquePostdeviant_ElectrodeChannel_08.png

2 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O25C1.mcd/Raster_NegativeSpikes_Triggers_UniquePostdeviant_ElectrodeChannel_09.png

1 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE85_31.07.2019_O25C1.mcd/Raster_NegativeSpikes_Triggers_Deviant_ElectrodeChannel_31.png

1 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE85_31.07.2019_O25C1.mcd/Raster_NegativeSpikes_Triggers_Deviant_ElectrodeChannel_32.png

2 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE83_29.07.2019_O25C2.mcd/Raster_NegativeSpikes_Triggers_Standard_ElectrodeChannel_03.png

1 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE83_29.07.2019_O25C2.mcd/Raster_NegativeSpikes_Triggers_Standard_ElectrodeChannel_04.png

1 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE83_29.07.2019_O25C2.mcd/Raster_NegativeSpikes_Triggers_Standard_ElectrodeChannel_05.png

1 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE83_29.07.2019_O25C2.mcd/Raster_NegativeSpikes_Triggers_Standard_ElectrodeChannel_06.png

2 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE83_29.07.2019_O25C2.mcd/Raster_NegativeSpikes_Triggers_Standard_ElectrodeChannel_07.png

2 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE83_29.07.2019_O25C2.mcd/Raster_NegativeSpikes_Triggers_Standard_ElectrodeChannel_08.png

2 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE83_29.07.2019_O25C2.mcd/Raster_NegativeSpikes_Triggers_Standard_ElectrodeChannel_09.png

2 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE83_29.07.2019_O25C2.mcd/Raster_NegativeSpikes_Triggers_UniquePredeviant_ElectrodeChannel_04.png

1 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE83_29.07.2019_O25C2.mcd/Raster_NegativeSpikes_Triggers_UniquePredeviant_ElectrodeChannel_05.png

1 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE83_29.07.2019_O25C2.mcd/Raster_NegativeSpikes_Triggers_UniquePredeviant_ElectrodeChannel_06.png

1 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE83_29.07.2019_O25C2.mcd/Raster_NegativeSpikes_Triggers_UniquePredeviant_ElectrodeChannel_08.png

2 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE83_29.07.2019_O25C2.mcd/Raster_NegativeSpikes_Triggers_UniquePostdeviant_ElectrodeChannel_05.png

2 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE83_29.07.2019_O25C2.mcd/Raster_NegativeSpikes_Triggers_UniquePostdeviant_ElectrodeChannel_06.png

1 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O25C2.mcd/Raster_NegativeSpikes_Triggers_Standard_ElectrodeChannel_31.png

1 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O25C2.mcd/Raster_NegativeSpikes_Triggers_UniquePredeviant_ElectrodeChannel_31.png

2 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O25C2.mcd/Raster_NegativeSpikes_Triggers_Deviant_ElectrodeChannel_08.png

2 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O25C2.mcd/Raster_NegativeSpikes_Triggers_UniquePostdeviant_ElectrodeChannel_31.png

1 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O25UC1.mcd/Raster_NegativeSpikes_Triggers_Standard_ElectrodeChannel_06.png

1 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O25UC1.mcd/Raster_NegativeSpikes_Triggers_Standard_ElectrodeChannel_07.png

1 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O25UC1.mcd/Raster_NegativeSpikes_Triggers_Standard_ElectrodeChannel_08.png

1 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O25UC1.mcd/Raster_NegativeSpikes_Triggers_Standard_ElectrodeChannel_09.png

1 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O25UC1.mcd/Raster_NegativeSpikes_Triggers_Predeviant_ElectrodeChannel_06.png

1 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O25UC1.mcd/Raster_NegativeSpikes_Triggers_Predeviant_ElectrodeChannel_07.png

1 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O25UC1.mcd/Raster_NegativeSpikes_Triggers_Predeviant_ElectrodeChannel_08.png

1 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O25UC1.mcd/Raster_NegativeSpikes_Triggers_Predeviant_ElectrodeChannel_09.png

1 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O25UC1.mcd/Raster_NegativeSpikes_Triggers_Postdeviant_ElectrodeChannel_06.png

1 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O25UC1.mcd/Raster_NegativeSpikes_Triggers_Postdeviant_ElectrodeChannel_07.png

1 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O25UC1.mcd/Raster_NegativeSpikes_Triggers_Postdeviant_ElectrodeChannel_08.png

1 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O25UC1.mcd/Raster_NegativeSpikes_Triggers_Postdeviant_ElectrodeChannel_09.png

2 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE85_31.07.2019_O25UC1.mcd/Raster_NegativeSpikes_Triggers_Deviant_ElectrodeChannel_31.png

2 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE85_31.07.2019_O25UC1.mcd/Raster_NegativeSpikes_Triggers_Deviant_ElectrodeChannel_32.png

1 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O25UC2.mcd/Raster_NegativeSpikes_Triggers_Deviant_ElectrodeChannel_08.png

2 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_O25UC2.mcd/Raster_NegativeSpikes_Triggers_Deviant_ElectrodeChannel_09.png

2 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_MS.mcd/Raster_NegativeSpikes_Triggers_C1_ElectrodeChannel_06.png

2 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_MS.mcd/Raster_NegativeSpikes_Triggers_C1_ElectrodeChannel_07.png

2 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_MS.mcd/Raster_NegativeSpikes_Triggers_C1_ElectrodeChannel_08.png

2 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_MS.mcd/Raster_NegativeSpikes_Triggers_C1_ElectrodeChannel_09.png

1 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_DOC1.mcd/Raster_NegativeSpikes_Triggers_Standard_ElectrodeChannel_06.png

1 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_DOC1.mcd/Raster_NegativeSpikes_Triggers_Standard_ElectrodeChannel_07.png

1 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_DOC1.mcd/Raster_NegativeSpikes_Triggers_Standard_ElectrodeChannel_08.png

1 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_DOC1.mcd/Raster_NegativeSpikes_Triggers_Standard_ElectrodeChannel_09.png

1 gwenview /media/loaloa/Samsung_T5/output_lowthr/MUA_output_lowthr/mGE84_30.07.2019_DOC1.mcd/Raster_NegativeSpikes_Triggers_Standard_ElectrodeChannel_11.png

one_weights = [2 ,2 ,2 ,2 ,2 ,2 ,2 ,2 ,2 ,2 ,2 ,2 ,1 ,2 ,1 ,1 ,1 ,1 ,2 ,1 ,1 ,1 ,1 ,2 ,1 ,1 ,2 ,1 ,1 ,2 ,1 ,1 ,1 ,2 ,2 ,2 ,2 ,1 ,1 ,1 ,2 ,2 ,1 ,1 ,2 ,2 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,2 ,2 ,1 ,2 ,2 ,2 ,2 ,2 ,1 ,1 ,1 ,1 ,1]

"""


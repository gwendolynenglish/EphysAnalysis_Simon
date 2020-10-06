# coding: utf-8
import MUA_constants as const
from preprocessing import compress_CSVs, process_data
import plotting


import os
from MUA_constants import P


"""Process data """
# process_data(how={'nbatches': 3, 'batch':0})
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
from onset_offset_classif import onset_offset_response
from onset_offset_classif import onset_offset_labels
from onset_offset_classif import lapl_kernel_SVM
from onset_offset_classif import get_onset_offset_classification
from onset_offset_classif import onoff_heatmap
from onset_offset_classif import onoff_barplot

# onset_offset_response(plots_dest_dir_appdx='../plots/onset_offset_chnls', csv_dest_dir_appdx='../onset_offset_model')
# onset_offset_response(plots_dest_dir_appdx='../plots/onset_offset_regions', single_channels=False, csv_dest_dir_appdx=None)
# onset_offset_labels(dest_dir_appdx='../onset_offset_model')

training_data_dir = '/mnt/Samsung_T5/output_lowthr/onset_offset_model'
# lapl_kernel_SVM(training_data_dir=training_data_dir, dest_dir_appdx='../onset_offset_model', parameter_search=True)
# lapl_kernel_SVM(training_data_dir=training_data_dir, dest_dir_appdx='../onset_offset_model', plot_cv_result=True)
# lapl_kernel_SVM(training_data_dir=training_data_dir, analyize_confusions=True)u

training_data_chnl_map_file = '/mnt/Samsung_T5/output_lowthr/chnls_map.csv'
MUA_output_data_chnl_map_file = '/mnt/Samsung_T5/output_gwen_data/channel_mappings/S1Th_LayerAssignment_22.10.19.csv'

# training data
# train_data = get_onset_offset_classification(which_data='training', training_data_dir=training_data_dir, 
#                                              training_data_chnl_map_file=training_data_chnl_map_file, MUA_output_data_chnl_map_file=MUA_output_data_chnl_map_file)
# onoff_heatmap(train_data, dest_dir_appdx='../onset_offset_model/plots', fig_height=11)
# train_data_all = get_onset_offset_classification(which_data='training', training_data_dir=training_data_dir, keep_labels=[0,1,3],
#                                                  training_data_chnl_map_file=training_data_chnl_map_file, MUA_output_data_chnl_map_file=MUA_output_data_chnl_map_file)
# onoff_barplot(train_data_all, dest_dir_appdx='../onset_offset_model/plots')









# change MUA_constans ALL_MICE and MUA init P["outputPath"] to new data
# dest_dir_appdx = '../pred_gwendata_highthr'
# dest_dir_plot_appdx = '../pred_gwendata_highthr/gwendata_plots'
# os.makedirs(f'{P["outputPath"]}/{dest_dir_plot_appdx}', exist_ok=True)

# new_data = get_onset_offset_classification(which_data='MUA_output', training_data_dir=training_data_dir, dest_dir_appdx=dest_dir_appdx, cached_prediction=True,
#                                        training_data_chnl_map_file=training_data_chnl_map_file, MUA_output_data_chnl_map_file=MUA_output_data_chnl_map_file)
# onoff_heatmap(new_data, dest_dir_appdx=dest_dir_plot_appdx, fig_height=14)
# new_data_all = get_onset_offset_classification(which_data='MUA_output', training_data_dir=training_data_dir, dest_dir_appdx=dest_dir_appdx, keep_labels=[0,1,2,3],
#                                                  training_data_chnl_map_file=training_data_chnl_map_file, MUA_output_data_chnl_map_file=MUA_output_data_chnl_map_file)
# onoff_barplot(new_data_all, dest_dir_appdx=dest_dir_plot_appdx)

# both training (but highthreshold, trained on low threshold) and new data
# change MUA_constans to all mice!
dest_dir_appdx = '../pred_alldata_highthr'
dest_dir_plot_appdx = '../pred_alldata_highthr/alldata_plots'
os.makedirs(f'{P["outputPath"]}/{dest_dir_plot_appdx}', exist_ok=True)

MUA_output_data_chnl_map_file = '/mnt/Samsung_T5/output_all_data_highthr/S1th_LayerAssignment_22.10.19.csv'


both_data = get_onset_offset_classification(which_data='MUA_output', training_data_dir=training_data_dir, dest_dir_appdx=dest_dir_appdx, cached_prediction=True, 
                                       training_data_chnl_map_file=training_data_chnl_map_file, MUA_output_data_chnl_map_file=MUA_output_data_chnl_map_file)
onoff_heatmap(both_data, dest_dir_appdx=dest_dir_plot_appdx, fig_height=15)

both_data = get_onset_offset_classification(which_data='MUA_output', training_data_dir=training_data_dir, dest_dir_appdx=dest_dir_appdx, keep_labels=[0,1,2,3],
                                                 training_data_chnl_map_file=training_data_chnl_map_file, MUA_output_data_chnl_map_file=MUA_output_data_chnl_map_file)
onoff_barplot(both_data, dest_dir_appdx=dest_dir_plot_appdx)













# _200_Raster_NegativeSpikes_Triggers_Deviant_ElectrodeChannel_21.png
# _231_Raster_NegativeSpikes_Triggers_C1_ElectrodeChannel_21.png
# _342_Raster_NegativeSpikes_Triggers_Deviant_ElectrodeChannel_12.png

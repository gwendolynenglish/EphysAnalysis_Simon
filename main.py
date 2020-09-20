# coding: utf-8
import MUA_constants as const
from preprocessing import compress_CSVs, process_data
import plotting


"""Process data """
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
from onset_offset_classif import onset_offset_response, onset_offset_labels, lapl_kernel_SVM, get_onset_offset_classification, onoff_heatmap, onoff_barplot
# onset_offset_response(plots_dest_dir_appdx='../plots/onset_offset_chnls', csv_dest_dir_appdx='../onset_offset_model')
# onset_offset_response(plots_dest_dir_appdx='../plots/onset_offset_regions', single_channels=False, csv_dest_dir_appdx=None)
# onset_offset_labels(dest_dir_appdx='../onset_offset_model')

training_data_dir = '/mnt/Samsung_T5/output_lowthr/onset_offset_model'
# lapl_kernel_SVM(training_data_dir=training_data_dir, dest_dir_appdx='../onset_offset_model', parameter_search=True)
# lapl_kernel_SVM(training_data_dir=training_data_dir, dest_dir_appdx='../onset_offset_model', plot_cv_result=True)
# lapl_kernel_SVM(training_data_dir=training_data_dir, analyize_confusions=True)u

training_data_chnl_map_file = f'{const.P["outputPath"]}/../../output_lowthr/chnls_map.csv'
MUA_output_data_chnl_map_file = f'{const.P["outputPath"]}/../channel_mappings/S1Th_LayerAssignment_22.10.19.csv'

# training data
# train_data = get_onset_offset_classification(which_data='training', training_data_dir=training_data_dir, 
#                                              training_data_chnl_map_file=training_data_chnl_map_file, MUA_output_data_chnl_map_file=MUA_output_data_chnl_map_file)
# onoff_heatmap(train_data, dest_dir_appdx='../onset_offset_model/plots', fig_height=11)
# train_data_all = get_onset_offset_classification(which_data='training', training_data_dir=training_data_dir, keep_labels=[0,1,3],
#                                                  training_data_chnl_map_file=training_data_chnl_map_file, MUA_output_data_chnl_map_file=MUA_output_data_chnl_map_file)
# onoff_barplot(train_data_all, dest_dir_appdx='../onset_offset_model/plots')


# change MUA_constans ALL_MICE and MUA init P["outputPath"] to new data
dest_dir_appdx = '../pred_gwendata_lowthr/gwendata_plots'
new_data = get_onset_offset_classification(which_data='MUA_output', training_data_dir=training_data_dir, dest_dir_appdx='../pred_gwendata_lowthr', 
                                       training_data_chnl_map_file=training_data_chnl_map_file, MUA_output_data_chnl_map_file=MUA_output_data_chnl_map_file)
onoff_heatmap(new_data, dest_dir_appdx=dest_dir_appdx, fig_height=14)
new_data_all = get_onset_offset_classification(which_data='MUA_output', training_data_dir=training_data_dir, dest_dir_appdx='../pred_gwendata_lowthr', keep_labels=[0,1,2,3],
                                                 training_data_chnl_map_file=training_data_chnl_map_file, MUA_output_data_chnl_map_file=MUA_output_data_chnl_map_file)
onoff_barplot(new_data_all, dest_dir_appdx=dest_dir_appdx)


# both training and new data
dest_dir_appdx = '../pred_gwendata_lowthr/gwendata_plus_train_plots'
both_data = get_onset_offset_classification(which_data='both', training_data_dir=training_data_dir, dest_dir_appdx='../pred_gwendata_lowthr', 
                                       training_data_chnl_map_file=training_data_chnl_map_file, MUA_output_data_chnl_map_file=MUA_output_data_chnl_map_file)
onoff_heatmap(both_data, dest_dir_appdx=dest_dir_appdx, fig_height=15)

both_data = get_onset_offset_classification(which_data='both', training_data_dir=training_data_dir, dest_dir_appdx='../pred_gwendata_lowthr', keep_labels=[0,1,2,3],
                                                 training_data_chnl_map_file=training_data_chnl_map_file, MUA_output_data_chnl_map_file=MUA_output_data_chnl_map_file)
onoff_barplot(both_data, dest_dir_appdx=dest_dir_appdx)













# _200_Raster_NegativeSpikes_Triggers_Deviant_ElectrodeChannel_21.png
# _231_Raster_NegativeSpikes_Triggers_C1_ElectrodeChannel_21.png
# _342_Raster_NegativeSpikes_Triggers_Deviant_ElectrodeChannel_12.png

%% Prepare spike tirain data (for Figure 3,4,5,6,7,8)
% need to change the folder adress into the correct address
% please open the program (auc_folders.m / fof_folders.m / mpfc_folders)

makeSpikeTrain('auc_folders');
makeSpikeTrain('fof_folders');
makeSpikeTrain('mpfc_folders');

%% Prepare intermediate data (for Figure 3,4,5)
analysis_dir1 = eval('auc_folders');
analysis_dir2 = eval('fof_folders');
analysis_dir3 = eval('mpfc_folders');
analysis_dir = [analysis_dir1;analysis_dir2;analysis_dir3];

for i=1:length(analysis_dir)
    pathname = analysis_dir{i};
    Task_kaiseki_tokyo1_20220201_make_ave_trace(pathname)
    Task_kaiseki_tokyo1_20220912_sound_choice_separate2_detail(pathname)
    norm_Task_kaiseki_tokyo1_20220316_sound_choice_separate2(pathname)
    norm_Task_kaiseki_tokyo1_20220316_sound_choice_separate2_anova(pathname)
    norm_Task_kaiseki_tokyo1_20220316_sound_choice_separate2_detail(pathname)
end

%% Prepare GLM encoding data (for Figure 6)
% need to change the folder adress into the correct address
% please open the program (auc_folders.m / fof_folders.m / mpfc_folders)

Task_kaiseki_tokyo_encode_20220912_all_all('auc_folders');
Task_kaiseki_tokyo_encode_20220912_all_all('fof_folders');
Task_kaiseki_tokyo_encode_20220912_all_all('mpfc_folders');

%% Prepare intermediate data (for Figure 8)
analysis_dir1 = eval('auc_folders');
analysis_dir2 = eval('fof_folders');
analysis_dir3 = eval('mpfc_folders');
analysis_dir = [analysis_dir1;analysis_dir2;analysis_dir3];

for i=1:length(analysis_dir)
    pathname = analysis_dir{i};
    State_dynamics_20220917_QR_CV_depth_control(pathname)
end

%%
% Below section for visualization 
% The run time of almost all scripts were within 1 minute (in 10 minutes at maximum).
%%
%% Make Figure 1 (No need to execute the above scripts)
Figure1b_Dual_behave_analysis_block1_200706_plot5_precise2;
% need to select the appropriate file
% please open the program 

Figure1c_Dual_Full_psycho_analysis_230630_all_session;
% need to change the folder adress into the correct address
% please open the program 

Figure1d_Dual_Full_psycho_analysis_220818_all;
% need to change the folder adress into the correct address
% please open the program 

Figure1e_Dual_behave_analysis_block1_220818_lick_long4;
% need to change the folder adress into the correct address
% please open the program 

Figure1f_Dual_behave_analysis_block1_220818_short_long;
% need to change the folder adress into the correct address
% please open the program 

%% Make Figure 2 (No need to execute the above scripts)
Figure2b_Dual_behave_analysis_block1_200706_plot5_precise2
% need to select the appropriate file
% please open the program 

Figure2b_Dual_RL_model_block1_220818_plot_values2
% need to select the appropriate file
% please open the program 

Figure2c_Dual_RL_analysis_220818_para_compare
% need to change the folder adress into the correct address
% please open the program 

Figure2d_Dual_RL_analysis_220314_para_all_220909_2
% need to change the folder adress into the correct address
% please open the program 

Figure2ef_process_20220818_prior_process_value_all2
% need to change the folder adress into the correct address
% please open the program 

%% Make Figure 3
Figure3b_process_20220516_norm_ave_trace_depth_control('mpfc_folders')
Figure3b_process_20220912_simple_process3_raw_depth_control('mpfc_folders')

Figure3c_Task_kaiseki_tokyo1_20220210_sound_trace
% need to change the folder adress into the correct address
% please open the program 

Figure3d_process_20230701_simple_process3_raw_depth_control('mpfc_folders')

Figure3e_Task_kaiseki_tokyo1_20220215_prior_trace
% need to change the folder adress into the correct address
% please open the program 

%% Make Figure 4
Figure4a_process_20230711_simple_process_test4_depth_control('mpfc_folders')
Figure4b_process_20220516_simple_process_test4_depth_control('mpfc_folders')
Figure4c_process_20220516_simple_process_test5_depth_control('mpfc_folders')
Figure4d_process_20240102_simple_process_test4_depth_control('mpfc_folders')

%% Make Figure 5
Figure_b_process_20220516_norm_ave_trace_depth_control('auc_folders')
Figure_b_process_20220912_simple_process3_raw_depth_control('auc_folders')

Figure_c_Task_kaiseki_tokyo1_20220210_sound_trace
% need to change the folder adress into the correct address
% please open the program 

Figure_d_process_20230701_simple_process3_raw_depth_control('auc_folders')

Figure_e_Task_kaiseki_tokyo1_20220215_prior_trace
% need to change the folder adress into the correct address
% please open the program 

Figure_f_process_20230711_simple_process_test4_depth_control('auc_folders')
Figure_g_process_20220516_simple_process_test4_depth_control('auc_folders')
Figure_h_process_20220516_simple_process_test5_depth_control('auc_folders')
Figure_i_process_20240102_simple_process_test4_depth_control('auc_folders')

%% Make Figure 6
Figure6ab_Task_20231002_encoding_model_raw
% need to change the folder adress into the correct address
% please open the program 

Figure6b_Task_20231002_encoding_model_raw_WO_sound
% need to change the folder adress into the correct address
% please open the program 

Figure6cd_process_20231002_encoding_model_compare_depth_control

%% Make Figure 7
Figure7a1_Population_decoder_20220520_SLR_process_one_session
% need to change the folder adress into the correct address
% please open the program 

Figure7a2_Population_decoder_20220520_SLR_process_one_session
% need to change the folder adress into the correct address
% please open the program 

Figure7a3_Population_decoder_SLR_one_session_prior
% need to change the folder adress into the correct address
% please open the program 

Figure7b_process_20220520_0912_SLR_limit_neuron_depth_control('auc_folders')
Figure7b_process_20220520_0912_SLR_limit_neuron_depth_control('fof_folders')
Figure7b_process_20220520_0912_SLR_limit_neuron_depth_control('mpfc_folders')

Figure7c_process_20230920_0915_SLR_long_short_integrate('auc_folders')
Figure7c_process_20230920_0915_SLR_long_short_integrate('fof_folders')
Figure7c_process_20230920_0915_SLR_long_short_integrate('mpfc_folders')

Figure7d_process_20220520_1206_SLR_compare

Figure7e_process_20220520_0912_SLR_limit_neuron_depth_control('auc_folders')
Figure7e_process_20220520_0912_SLR_limit_neuron_depth_control('fof_folders')
Figure7e_process_20220520_0912_SLR_limit_neuron_depth_control('mpfc_folders')

Figure7e_process_20220520_0912_SLR_long_detail_limit_neuron2('auc_folders')
Figure7e_process_20220520_0912_SLR_long_detail_limit_neuron2('fof_folders')
Figure7e_process_20220520_0912_SLR_long_detail_limit_neuron2('mpfc_folders')

%% Make Figure 8
Figure8cde_process_20220917_state_dynamics_CV_depth_control('auc_folders')
Figure8cde_process_20220917_state_dynamics_CV_depth_control('fof_folders')
Figure8cde_process_20220917_state_dynamics_CV_depth_control('mpfc_folders')

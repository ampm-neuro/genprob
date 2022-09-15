

% for time warping
tses = [0.2 1.0 1.0 2.0 3.0 2.0];

% warp trial traces
[warp_traces, warp_frame_times, warp_trl_idx, warp_trl_mtx] = ...
    timewarp_traces(traces, frame_times, trl_idx, trl_mtx, tses);

% plot activity on every trial, one fig per neuron
tw_activity_plot_trial_hm(trl_mtx, trl_idx, medass_cell, ...
    frame_times, traces, 1:size(traces,1), 1:size(trl_mtx,1), 4, [5 25]);
tw_activity_plot_trial_hm_warp(trl_mtx, trl_idx, frame_times, traces, 1:size(traces,1), 1:size(trl_mtx,1), ceil([sum(tses(1:3)) sum(tses(4:end))]), tses);

% activity matrix (trials x time x neuron) for single session
trial_activity_mtx = image_trial_activity_mtx(session_path);

% plot with rows sorted by tone
all_matrices = tw_activity_plot_tone_hm_warp(trl_mtx, trl_idx, medass_cell, frame_times, traces,1:size(traces,1), 1:size(trl_mtx,1), 4, [5 25], tses);

% plot sorted cell means over warped trial
[~, srt_idx] = image_mean_activity_timewarp(...
    trl_mtx, trl_idx, frame_times, traces, 1:size(traces,1), unique(trl_idx), tses);

% make gif of single cell all trial plots
makeGif

% compute summary of correlation between activity on probe trials and how
% long the animal waited
ALL_image_mean_activity_timewarp_corr({'postprobe'}, 'all');

% syncronous activity
image_plot_trialtraces_trl_sync(trl_mtx, frame_times, trl_idx, traces, 1:size(traces,1), 30, srt_idx);

% cell reg nested functions
origin_fp_corrected = footprints_meta_reverse(cregstruct_meta, cregstruct_intermediates);
comp_footprints(old_fp, new_fp); % plots

% prep footprints
save_footprints_cellregprep;



%% plots of footprints
%load('C:\Users\ampm1\Desktop\cellreg_temp\658648m2\save_folder_all\aligned_data_struct.mat') % aligned_data_struct
trace_footprints(footprint_mtx) %single session)
trace_footprints_sessions % two sessions overlaid
trace_footprints_sessions_overlap % two sessions overlaid with distinct color for shared cells
footprints_cell = aligned_data_struct.spatial_footprints_corrected;
footprints_cell_chron = footprints_cell(session_chron_reorder_idx);
figure; trace_footprints_sessions_overlap(aligned_data_struct.spatial_footprints_corrected([2 3]), cell_regist_mtx(:,[2 3]))
all_footprints_colored(aligned_data_struct.spatial_footprints_corrected, cell_regist_mtx)


% plot and color footprints by activity peak time (uninteresting)
session_chron_reorder_idx = [1 9 15 2 10 16 3 11 17 4 12 18 5 13 19 6 14 20 7 8];
footprint_mtx = aligned_data_struct.spatial_footprints_corrected(session_chron_reorder_idx);
trace_footprint_activityPeak(aligned_data_struct.spatial_footprints_corrected(session_num), trl_mtx, trl_idx, frame_times, traces, tses)

save_footprints_cellregprep(mouse, folder)

% smooth footprints
all_footprints = get_file_paths_targeted('C:\Users\ampm1\Desktop\cellreg_temp\spatial footprints', 'footprints', 'prob');
for ifp = 15:size(all_footprints,1)
    load(all_footprints{ifp}, 'footprints');
    for ineuron = 1:size(footprints,1)
        footprints(ineuron, :,:) = smooth2a(squeeze(footprints(ineuron,:,:)),3);
    end
    save(all_footprints{ifp}, 'footprints');
end


% all footprint comparisons (single subject)
all_footprints_bottomUP_colored(aligned_data_struct.spatial_footprints_corrected, cell_regist_mtx)
all_footprints_topDown_colored(aligned_data_struct.spatial_footprints_corrected, cell_regist_mtx)

% proportion of cells in each session that were not active in any
% prior session
[props_vect] = prop_new(crm);
figure; errorbar_plot(mat2cell([prop_new(cell_regist_mtx_658648m2);[prop_new(cell_regist_mtx_651049m1)]], 2, ones(1,20) ));


%% matrix of overlapping cells
%
mevar_subjs = [{'651049m1'} {'658648m2'} {'691359m1'} {'690330m1'} {'691359m2'} {'696944m3'} {'699437m4'} {'699438m3'}];
hivar_subjs = [{'687034m4'} {'683472m2'} {'699437m2'} {'699438m1'} {'683472m3'}];

overlap_mtx_mevar = overlapping_population_average(mevar_subjs);
overlap_mtx_hivar = overlapping_population_average(hivar_subjs);
%{
cell_regist_mtx = cell_registered_struct.cell_to_index_map;
session_chron_reorder_idx = [1 9 15 2 10 16 3 11 17 4 12 18 5 13 19 6 14 20 7 8];
cell_overlap_mtx({cell_regist_mtx}, session_chron_reorder_idx);
%}

% plot of overlap vs probe similarity
pairwise_probesimilarity_vs_overlap


%% flow of cells between active and inactive over imaging sessions
sankey_local(cell_regist_mtx, 1:8)
sankey_full(cell_regist_mtx, 1:8)
sankey_full_probeColor(cell_regist_mtx, 1:8)

%% temporal correlations between common cells

% Compare stability of activity patterns over all days, subj level
all_subjects_cellcorr 

% Compare stability patterns over select days, cell level
pairwise_session_remapping_wrap 

% correlate between neuronal stability and probe stability 
cellCorrs_v_probeCorrs



%% mean normalized rate of all neurons in each stage
[mean_rate_matrix] = mean_rate_vect('mevar', [{'651049m1','658648m2', '690330m1', '690330m2', '691359m1'}]);
figure; imagesc(mean_rate_matrix)


%% mds summary



%% distances
[activity_mtx, session_mtx_cell, trial_mtx_cell_means, session_mtx_cell_means, session_number_idx, trial_number_idx, time_bin_idx] = image_subject_allsesh_decode('683472m2', cell_regist_mtx, 1:8, tses);%distances = distance_window(session_mtx_cell_means', session_mtx_cell_means', 2);
image_subject_allsesh_decode('683472m2', cell_regist_mtx, 1:8, tses);
[mds_coords, stress, distances_mtx] = mds_plot(activity_mtx, 2, session_number_idx);
% [mds_coords, stress, distances_mtx] = mds_plot(session_mtx_cell_means, 2, session_number_idx);

% moment to moment distances
% ouput from: [activity_mtx, session_mtx_cell, trial_mtx_cell_means, session_mtx_cell_means, session_number_idx, trial_number_idx, time_bin_idx] = image_subject_allsesh_decode('658648m2', cell_regist_mtx, 1:8, tses);
current_session = 3;
comparison_sessions = [2 4];
for itrl = 1:2
    [all_distances, mean_diff_distance] = ...
        image_timeBin_popDist_line(current_session, comparison_sessions, itrl, activity_mtx, cell_regist_mtx, session_number_idx, trial_number_idx, time_bin_idx, tses); 
end

% correlate probe wait times with difference between population distances
% ouput from: [activity_mtx, session_mtx_cell, trial_mtx_cell_means, session_mtx_cell_means, session_number_idx, trial_number_idx, time_bin_idx] = image_subject_allsesh_decode('658648m2', cell_regist_mtx, 1:8, tses);
image_corr_popDist_waits('658648m2', 6, [4 7], activity_mtx, cell_regist_mtx, session_number_idx, trial_number_idx, time_bin_idx, tses)

% correlate as above but with r values instead of distances
image_corr_popCorr_waits('658648m2', 5, [4 3], cell_regist_mtx, tses)

% similarity to prior probe
image_similarity_to_prior_probe('658648m2', cell_regist_mtx, tses)



%% tracking cell reactivation over all 20 sessions
session_chron_reorder_idx = [1 9 15 2 10 16 3 11 17 4 12 18 5 13 19 6 14 20 7 8];
super_crm_out = super_crm({'651049m1','658648m2', '690330m1', '690330m2', '691359m1'});
image_engram_integration(super_crm_out(:, session_chron_reorder_idx))
%image_engram_integration_multi(crm_cell);

% plot all trials (across sessions) from single neurons
mouse = '683472m3'; 

load(['cell_reg_' mouse '.mat']); cell_regist_mtx = cell_registered_struct.cell_to_index_map; % cell_regist_mtx

last_prob_sessions = [];
for iprob = 1:6
prob_sessions = get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_imaging_hpc\' mouse ], ['var0' num2str(iprob)], 'LED');
prob_sessions = get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_hivar_imaging_hpc\' mouse ], ['var0' num2str(iprob)], 'LED');

last_prob_sessions = [last_prob_sessions; prob_sessions(end)];

end

session_cell = [...
get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_imaging_hpc\' mouse ], 'preprobe', 'LED');...
get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_imaging_hpc\' mouse ] , 'postprobe', 'LED');...
get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_imaging_hpc\' mouse ], 'var0', '-01', 'LED');...
get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_hivar_imaging_hpc\' mouse ], 'preprobe', 'LED');...
get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_hivar_imaging_hpc\' mouse ] , 'postprobe', 'LED');...
get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_hivar_imaging_hpc\' mouse ], 'var0', '-01', 'LED');...
last_prob_sessions];

session_chron_reorder_idx = [1 9 15 2 10 16 3 11 17 4 12 18 5 13 19 6 14 20 7 8];
tses = [0.2 1.1 1.0 2.0 8.5 2.0];
neuron = 15;
image_warped_singlecell_trials_hm_multisesh(neuron, cell_regist_mtx, session_cell, session_chron_reorder_idx, tses)

    % combine into one figure
    load('image_pop_mds_full_save.mat', 'session_mtx_cell')
    % or
    [activity_mtx, session_mtx_cell, trial_mtx_cell_means, session_mtx_cell_means, session_number_idx, trial_number_idx, time_bin_idx] = image_subject_allsesh_decode(mouse, cell_regist_mtx, 1:20, tses)
    
    min_active_sessions = 17;
    cells_with_min_active_sessions = find(sum(cell_regist_mtx>0,2)>=min_active_sessions);
    for ic = cells_with_min_active_sessions'
        all_trials_one_cell_fr(ic, session_mtx_cell, tses);
        title(num2str(ic))
    end
    
    % add footprint tracking
    cell_num = 197;
    figure; sessions = [14:20] ; 
    for isesh = 1:length(sessions) 
        disp(num2str(cell_regist_mtx(cell_num, sessions(isesh)))); 
        subplot(1,length(sessions),isesh); 
        trace_footprints(footprint_mtx{sessions(isesh)}); 
        hold on; 
        trace_footprint(squeeze(footprint_mtx{sessions(isesh)}(cell_regist_mtx(cell_num, sessions(isesh)),:,:)), [0 0 0]);
    end
    

% firing rate field metrics
all_info_content = all_cells_infoContent(1:646, session_mtx_cell); % dispersion, not info content
reliability_corm = ALL_reliability_corm(session_mtx_cell, cell_regist_mtx); % plot average cell matrix of trial by trial correlation of activity

% comparing info content of top down , bottom up , and misc cells
[topDown_idx, bottomUp_idx] = cell_type_idx(cell_regist_mtx(:,session_chron_reorder_idx));
[supervect, supervect_idx] = image_compare_cellType_activity(session_cell(session_chron_reorder_idx), cell_regist_mtx(:,session_chron_reorder_idx), tses);
[all_specificity_mevar, all_reliability_mevar, all_reactivationCt_mevar, errorbarplot_in_mevar] = image_compare_cellType_activity_multi('mevar', [{'651049m1'} {'658648m2'}]); 
[all_specificity_hivar, all_reliability_hivar, all_reactivationCt_hivar, errorbarplot_in_hivar] = image_compare_cellType_activity_multi('hivar', [{'683472m2'} {'683472m3'}]);
    % can combine output from mevar and hivar to plot both^^


    
%% relationship between population activity, learning, and integration
cellCorrs_v_probeCorrs_priorprobe_v_lastproblem
% also see: for behavior only correlations (can be used with imaging mice
% or behavior mice)
cellCorrs_v_probeCorrs_priorprobe_v_lastproblem_behonly


% pca attempt
%image_pca_plot_session(all_matrices, tses)
image_pca_plot_session2 % do this but with problem days
image_pca_plot_session_trialorder


% prediction scores vs overlap !!!
predict_v_overlap




% desired subjects
mouse_ids = [{'658648m2'} {'651049m1'}];
%mouse_ids = [{'683472m2'} {'683472m3'}];

% misc
tses = [0.2 1.1 1.0 2.0 8.5 2.0];

% for each subj
all_activity_mtx = [];
for isubj = 1:length(mouse_ids)
isubj
    
    % load cell_regist_mtx
    load(['cell_reg_' mouse_ids{isubj} '.mat'])
    cell_regist_mtx = cell_registered_struct.cell_to_index_map; 
    
    % compute activity matrix
    %[activity_mtx, session_mtx_cell, session_number_idx, trial_number_idx, time_bin_idx] = ...
    %    image_subject_allsesh_decode_bytrial(mouse_ids{isubj}, cell_regist_mtx, 1:20, tses);
    [activity_mtx, session_mtx_cell, session_number_idx, trial_number_idx, time_bin_idx] = ...
        image_subject_allsesh_decode_bytrial(mouse_ids{isubj}, cell_regist_mtx, [1:3:19 20], tses);

    % trials*sessions, time*neuron
    all_activity_mtx = [all_activity_mtx activity_mtx];
    
end

save('mds_input.mat', '-v7.3')


%% combine matrices
load('mds_input.mat')
downsamp = 1;
[mds_coords, stress, distances_mtx] = mds_plot(all_activity_mtx, 2, session_number_idx, downsamp);
[mds_coords, stress] = mds_plot_2(distances_mtx, 2, session_number_idx, downsamp);
function [activity_mtx, activity_mtx_pca, hm_cells, session_number_idx, trial_number_idx, time_bin_idx] = image_subject_allsesh_PCA(subject_id, cell_reg_mtx, sessions, time_series_event_spacing)
% perform PCA on activity of all neurons over all trials from all sessions
% plot first 3 components



%% get sessions for subject
session_path = ['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_imaging_hpc\' subject_id];
all_preprobes = get_file_paths_targeted(session_path, 'preprobe', 'LED');
    preprobe_sessions = sessions(sessions<7);
    kept_preprobes = [];
    for is = preprobe_sessions
        kept_preprobes = [kept_preprobes; all_preprobes(contains(all_preprobes, ['probe_0' num2str(is)]))];
    end
all_postprobes = get_file_paths_targeted(session_path, 'postprobe', 'LED');
    postprobe_sessions = sessions(sessions>=7)-6;
    kept_postprobes = [];
    for is = postprobe_sessions
        kept_postprobes = [kept_postprobes; all_postprobes(contains(all_postprobes, ['probe_0' num2str(is)]))];
    end
session_list = [kept_preprobes; kept_postprobes];



%% Compute and activity matrices

% preallocate
session_mtx_cell = cell(size(session_list,1), 1);

% iterate through sessions
for isesh = 1:size(session_list,1)
    
    % load session
    load(session_list{isesh})

    % only reward-unavailable trials
    hm_cell_trials = unique(trl_idx);
    hm_cell_trials = intersect(hm_cell_trials, find(trl_mtx(:,3)==0));
    
    % session matrix for every cell
    hm_cells = tw_activity_trial_hm_full_warp(trl_mtx, trl_idx, frame_times, traces, 1:size(traces,1), hm_cell_trials, [3 12], time_series_event_spacing);
    
    % remove times bins with nans on any trial
    for icell = 1:length(hm_cells)
        hm_cells{icell} = hm_cells{icell}(:,~isnan(sum(hm_cells{icell},1)));
        %hm_cells{icell} = reshape(zscore_mtx(hm_cells{icell}(:)), size(hm_cells{icell}));
    end

    % load
    session_mtx_cell{isesh} = hm_cells;
    
end



%% Reorganize sessions to match universal cell registry

% iterate through sessions
for isesh = 1:size(session_list,1)
    
    % matrix of zeros for inactive cells
    inactive_mtx = zeros(size(session_mtx_cell{isesh}{1}));

    % session matrix with all cells
    session_mtx_cell_hold = cell(size(cell_reg_mtx,1),1);
    
    % preload all cells as inactive
    for icell = 1:length(cell_reg_mtx(:,sessions(isesh)))
        session_mtx_cell_hold{icell} = inactive_mtx;
    end
    
    % load active cells
    for icell = unique(cell_reg_mtx(cell_reg_mtx(:,sessions(isesh))>0,sessions(isesh)))'
    	session_mtx_cell_hold{find(cell_reg_mtx(:,sessions(isesh))==icell, 1)} = session_mtx_cell{isesh}{icell};
    end
    
    % load updated session matrix
    session_mtx_cell{isesh} = session_mtx_cell_hold;
    
end



%% Average activity on each trial 
%
% iterate through sessions
for isesh = 1:length(session_mtx_cell)
    
    % iterate through cells
    for icell = 1:length(session_mtx_cell{isesh})
        
        % average activity on each trial
        session_mtx_cell{isesh}{icell} = mean(session_mtx_cell{isesh}{icell},2);
        
    end
end
%}



%% Form indices for plotting

% preallocate session-level indices
session_number_idx = cell(size(session_mtx_cell));
trial_number_idx = cell(size(session_mtx_cell));
time_bin_idx = cell(size(session_mtx_cell));

% iterate through all sessions
for isesh = 1:length(session_mtx_cell)
    
    % preallocate cell-level indices
    session_number_idx{isesh} = cell(size(session_mtx_cell{isesh}));
    trial_number_idx{isesh} = cell(size(session_mtx_cell{isesh}));
    time_bin_idx{isesh} = cell(size(session_mtx_cell{isesh}));
    
    % iterate through all cells
    for icell = 1:length(session_mtx_cell{isesh})
        
        % key values
        num_trials = size(session_mtx_cell{isesh}{icell},1);
        num_time_bins = size(session_mtx_cell{isesh}{icell},2);
        
        % load session index
        session_number_idx{isesh}{icell} = repmat(isesh, size(session_mtx_cell{isesh}{icell}));
        
        % load trial index
        trial_number_idx{isesh}{icell} = repmat((1:num_trials)', 1, num_time_bins);
        
        % load time index
        time_bin_idx{isesh}{icell} = repmat(1:num_time_bins, num_trials, 1);
        
    end
end



%% Merge activity matrices

% merge neural activity
activity_mtx = cells_to_mtx(session_mtx_cell);

% merge session index
session_number_idx = cells_to_mtx(session_number_idx);
session_number_idx = session_number_idx(:,1);

% merge trial index
trial_number_idx = cells_to_mtx(trial_number_idx);
trial_number_idx = trial_number_idx(:,1);

% merge time index
time_bin_idx = cells_to_mtx(time_bin_idx);
time_bin_idx = time_bin_idx(:,1);



%% compute dimensionality reduction
[activity_mtx_pca, ~, ~, ~, explained] = pca(activity_mtx');
%nnmfactors = nnmf(activity_mtx, 3);

%% plot

% colors
colors = winter(length(session_mtx_cell));

% pcas 4 5 6
figure; hold on
pca_plot3(session_mtx_cell, activity_mtx_pca, session_number_idx, sessions, colors, [4 5 6])

% pcas 1 2 3
figure; hold on
pca_plot3(session_mtx_cell, activity_mtx_pca, session_number_idx, sessions, colors, [1 2 3])

%{
% pcas 3 4
figure; hold on
pca_plot2(session_mtx_cell, activity_mtx_pca, session_number_idx, sessions, colors, [3 4])

% pcas 2 3
figure; hold on
pca_plot2(session_mtx_cell, activity_mtx_pca, session_number_idx, sessions, colors, [2 3])

% pcas 1 2
figure; hold on
pca_plot2(session_mtx_cell, activity_mtx_pca, session_number_idx, sessions, colors, [1 2])


% nnmfs 4 5 6
figure; hold on
pca_plot3(session_mtx_cell, activity_mtx_pca, session_number_idx, sessions, colors, [4 5 6])


% nnmfactors 1 2 3
figure; hold on
pca_plot3(session_mtx_cell, nnmfactors, session_number_idx, sessions, colors, [1 2 3])
title('nnmfactors1')
%}


end



%% internal functions

function mtx_out = cells_to_mtx(nested_cells)
% function for merging the nested cells into a single matrix 
% (time bins, neurons; see also indices)

% key values
num_cells = length(nested_cells{1});
num_sessions = length(nested_cells);
num_timebins = size(nested_cells{1}{1},2);
num_trials_fn = nan(num_sessions,1);
for isesh_fn = 1:num_sessions
    num_trials_fn(isesh_fn) = size(nested_cells{isesh_fn}{1},1);
end
bins_per_sesh = num_timebins.*num_trials_fn;

% preallocate
mtx_out = nan(sum(bins_per_sesh), num_cells);


% iterate through all sessions
lo_idx = 1;
for isesh_fn = 1:length(nested_cells)

    % iterate through all cells
    for icell_fn = 1:length(nested_cells{isesh_fn})

        % rotate
        current_cell_activity = nested_cells{isesh_fn}{icell_fn}';
        
        % load
        mtx_out(lo_idx:lo_idx+bins_per_sesh(isesh_fn)-1,icell_fn) = current_cell_activity(:);

    end
    
    % update idx
    lo_idx = lo_idx+bins_per_sesh(isesh_fn);
end
end

function pca_plot3(session_mtx, pca_mtx, session_idx, session_input, colors, pc3)

    % sesh time bins and means
    session_means = nan(length(session_mtx), size(pca_mtx,2));
    for isesh = 1:length(session_mtx)

        % plot trial time bins
        sesh_idx = session_idx==isesh;
        plot3(pca_mtx(sesh_idx,pc3(1)), pca_mtx(sesh_idx,pc3(2)), pca_mtx(sesh_idx,pc3(3)), '.', 'color', colors(isesh,:))

        % load session means
        session_means(isesh,:) = mean(pca_mtx(sesh_idx,:));
    end

    % plot line between session means
    plot3(session_means(:,pc3(1)), session_means(:,pc3(2)), session_means(:,pc3(3)), 'k-', 'linewidth', 3)

    % aesthetics
    xlabel(['PC ' num2str(pc3(1))])
    ylabel(['PC ' num2str(pc3(2))])
    zlabel(['PC ' num2str(pc3(3))])

    % legend
    legend_strings = [];
    for isesh = 1:length(session_mtx)
       legend_strings = [legend_strings; {['Probe ' num2str(session_input(isesh))]}];
    end
    legend(legend_strings)

end


function pca_plot2(session_mtx, pca_mtx, session_idx, session_input, colors, pc2)

    % sesh time bins and means
    session_means = nan(length(session_mtx), size(pca_mtx,2));
    for isesh = 1:length(session_mtx)

        % plot trial time bins
        sesh_idx = session_idx==isesh;
        plot(pca_mtx(sesh_idx,pc2(1)), pca_mtx(sesh_idx,pc2(2)), '.', 'color', colors(isesh,:))

        % load session means
        session_means(isesh,:) = mean(pca_mtx(sesh_idx,:));
    end

    % plot line between session means
    plot(session_means(:,pc2(1)), session_means(:,pc2(2)), 'k-', 'linewidth', 3)

    % aesthetics
    xlabel(['PC ' num2str(pc2(1))])
    ylabel(['PC ' num2str(pc2(2))])

    % legend
    legend_strings = [];
    for isesh = 1:length(session_mtx)
       legend_strings = [legend_strings; {['Probe ' num2str(session_input(isesh))]}];
    end
    legend(legend_strings)

end







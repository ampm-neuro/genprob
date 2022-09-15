function cell_turnover_timewarp_trials_halves(session_cell, cell_registration, time_series_event_spacing)
% plots average timewarped trials for every registered cell in every
% session, with inactive cells receiving zeros for activity



%% compute session matrix for each session
session_mtx_cell = cell(2, size(session_cell, 1));
session_max_cell = cell(2, size(session_cell, 1));
for isesh = 1:length(session_cell)
   
    % load session
    load(session_cell{isesh})
    
    % trials in first and second half
    all_trials = unique(trl_idx);
    first_half_trials = all_trials(all_trials<median(all_trials));
    second_half_trials = all_trials(all_trials>=median(all_trials));
    
    % session matrices
    [sorted_session_mtx_1sthalf, srt_idx_1sthalf, ~, max_time_bin_unsorted] = image_mean_activity_timewarp(...
    trl_mtx, trl_idx, frame_times, traces, 1:size(traces,1), first_half_trials, time_series_event_spacing);

    [sorted_session_mtx_2ndhalf, srt_idx_2ndhalf, ~, max_time_bin_unsorted] = image_mean_activity_timewarp(...
    trl_mtx, trl_idx, frame_times, traces, 1:size(traces,1), second_half_trials, time_series_event_spacing);

    % unsort
    unsort_idx = nan(size(srt_idx_1sthalf));
    unsort_idx(srt_idx_1sthalf) = 1:length(srt_idx_1sthalf);
    session_mtx_cell{1, isesh} = sorted_session_mtx_1sthalf(unsort_idx,:);
    
    unsort_idx = nan(size(srt_idx_2ndhalf));
    unsort_idx(srt_idx_2ndhalf) = 1:length(srt_idx_2ndhalf);
    session_mtx_cell{2, isesh} = sorted_session_mtx_2ndhalf(unsort_idx,:);
    
end


%% full session matrices
full_session_mtx_3d_1sthalf = nan(size(cell_registration,1), size(session_mtx_cell{isesh},2), length(session_mtx_cell));
full_session_mtx_3d_2ndhalf = nan(size(cell_registration,1), size(session_mtx_cell{isesh},2), length(session_mtx_cell));
for isesh = 1:length(session_mtx_cell)
    
    % preallocate halves
    full_session_mtx_1sthalf = zeros(size(cell_registration,1), size(session_mtx_cell{isesh},2));
    full_session_mtx_2ndhalf = zeros(size(cell_registration,1), size(session_mtx_cell{isesh},2));
    
    % indices
    active_cells_idx = cell_registration(:,isesh)>0;
    active_cell_session_ids = cell_registration(active_cells_idx,isesh);
    
    % load
    full_session_mtx_1sthalf(active_cells_idx, :) = session_mtx_cell{1, isesh}(active_cell_session_ids, :);
    full_session_mtx_3d_1sthalf(:,:,isesh) = full_session_mtx_1sthalf;
    
    full_session_mtx_2ndhalf(active_cells_idx, :) = session_mtx_cell{2, isesh}(active_cell_session_ids, :);
    full_session_mtx_3d_2ndhalf(:,:,isesh) = full_session_mtx_2ndhalf;
    
    
end



%% correlation between sessions using mutually active cells
% all cells is identical because correlations with a vector of zeros are
% always NaN

% preallocate
cell_corrs_mutual = cell(size(full_session_mtx_3d_1sthalf,3), size(full_session_mtx_3d_1sthalf,3));
cell_corrs_mutual_means = nan(size(full_session_mtx_3d_1sthalf,3), size(full_session_mtx_3d_1sthalf,3));

% iterate through session comparisons
for isesh01 = 1:size(full_session_mtx_3d_1sthalf,3)
    for isesh02 = 1:size(full_session_mtx_3d_2ndhalf,3)
        
        % active session comparison
        session_mtx01 = full_session_mtx_3d_1sthalf(:,:, isesh01);
        session_mtx02 = full_session_mtx_3d_2ndhalf(:,:, isesh02);
        
        % mutual cell idx
        mc_idx = sum(session_mtx01,2)>0 & sum(session_mtx01,2)>0;
        
        % mutual cells
        mc_session01 = session_mtx01(mc_idx,:);
        mc_session02 = session_mtx02(mc_idx,:);
        
        % iterate through cells
        cell_corrs = nan(size(mc_session01,1), 1);
        for icell = 1:size(mc_session01, 1)
            cell_corrs(icell) = corr(mc_session01(icell,:)', mc_session02(icell,:)');
        end
        
        % load overall
        cell_corrs_mutual{isesh01, isesh02} = cell_corrs;
        cell_corrs_mutual_means(isesh01, isesh02) = nanmean(cell_corrs);
    end
end

% plot means
figure;
ampm_pcolor(cell_corrs_mutual_means)
caxis([-1 1])
colorbar
xlabel('Probe first half'); ylabel('Probe second half')
title(remove_underscore('cell_corrs_mutual_means'))



%% informaiton content change of individual cells that are held across
% sessions






%% average activity of cells that are held across sessions?








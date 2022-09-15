function [activity_mtx, session_mtx_cell, trial_mtx_cell_means, session_mtx_cell_means, session_number_idx, trial_number_idx, time_bin_idx] = image_subject_allsesh_decode(session_cell, cell_reg_mtx, time_series_event_spacing)
% perform PCA on activity of all neurons over all trials from all sessions
% plot first 3 components



%% Compute and activity matrices

% preallocate
session_mtx_cell = cell(size(session_cell,1), 1);

% iterate through sessions
for isesh = 1:size(session_cell,1)
    
    % load session
    load(session_cell{isesh})
    
    hm_cell_trials = unique(trl_idx);
    %hm_cell_trials = intersect(hm_cell_trials, find(trl_mtx(:,3)==0)); % only reward-unavailable trials
    %hm_cell_trials = intersect(hm_cell_trials, find(trl_mtx(:,3)==0)); % only reward-available trials
    
    % session matrix for every cell
    hm_cells = tw_activity_trial_hm_full_warp(trl_mtx, trl_idx, frame_times, traces, 1:size(traces,1), hm_cell_trials, [3 12], time_series_event_spacing);
    
    % remove times bins with nans on every trial
    for icell = 1:length(hm_cells)
        hm_cells{icell} = hm_cells{icell}(:,sum(~isnan(hm_cells{icell}))>0);
    end
    %}

    % load
    session_mtx_cell{isesh} = hm_cells;
end

smc_length = [];
for isesh = 1:length(session_mtx_cell)
    smc_length = [smc_length length(session_mtx_cell{isesh})];
end
%smc_length = smc_length;
%crm_length = sum(cell_reg_mtx>0);

%% Reorganize sessions to match universal cell registry

% iterate through sessions
for isesh = 1:size(session_cell,1)

    % matrix of zeros for inactive cells
    inactive_mtx = zeros(size(session_mtx_cell{isesh}{1}));

    % session matrix with all cells
    session_mtx_cell_hold = cell(size(cell_reg_mtx,1),1);
    
    % preload all cells as inactive
    for icell = 1:length(cell_reg_mtx(:,isesh))
        session_mtx_cell_hold{icell} = inactive_mtx;
    end
    
    % load active cells
    for icell = unique(cell_reg_mtx(cell_reg_mtx(:,isesh)>0, isesh))'
    	session_mtx_cell_hold{find(cell_reg_mtx(:,isesh)==icell, 1)} = session_mtx_cell{isesh}{icell};
    end
    
    % load updated session matrix
    session_mtx_cell{isesh} = session_mtx_cell_hold;
    
end



%% Average activity 

% average session activity
session_mtx_cell_means = session_mtx_cell;
% iterate through sessions
for isesh = 1:length(session_mtx_cell_means)
    
    % iterate through cells
    for icell = 1:length(session_mtx_cell_means{isesh})
        
        % average activity on each trial
        session_mtx_cell_means{isesh}{icell} = mean(session_mtx_cell_means{isesh}{icell},2);
        
    end
    
    % average activity from each session
    session_mtx_cell_means{isesh} = mean(cell2mat(session_mtx_cell_means{isesh}'));
end

% matrix: sessions x cells
session_mtx_cell_means = cell2mat(session_mtx_cell_means);


% average trial activity
trial_mtx_cell_means = session_mtx_cell;
% iterate through sessions
for isesh = 1:length(trial_mtx_cell_means)
    
    % iterate through cells
    for icell = 1:length(trial_mtx_cell_means{isesh})
        
        % average activity on each trial
        trial_mtx_cell_means{isesh}{icell} = mean(trial_mtx_cell_means{isesh}{icell},2);
        
    end
    
    % average activity from each session
    trial_mtx_cell_means{isesh} = cell2mat(trial_mtx_cell_means{isesh}');
end

% matrix: trials x cells
trial_mtx_cell_means = cell2mat(trial_mtx_cell_means);




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

function internal_plot(session_mtx, pca_mtx, session_idx, session_input, colors, pc3)

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










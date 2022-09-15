function [activity_mtx, session_mtx_cell, session_number_idx, trial_number_idx, time_bin_idx] = image_subject_allsesh_decode_bytrial(subject_id, cell_reg_mtx, sessions, time_series_event_spacing)
% output activity_mtx is trials*sessions, time*neuron



%% get sessions for subject
%
% hivar v mevar
if ~isempty(get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_imaging_hpc\' subject_id], 'preprobe', 'LED'))
    session_path = ['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_imaging_hpc\' subject_id];
else
    session_path = ['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_hivar_imaging_hpc\' subject_id];
end



% all sessions
all_preprobes = get_file_paths_targeted(session_path, 'preprobe', 'LED');
    preprobe_sessions = 1:6;
    kept_preprobes = [];
    for is = preprobe_sessions
        kept_preprobes = [kept_preprobes; all_preprobes(contains(all_preprobes, ['probe_0' num2str(is)]))];
    end
all_postprobes = get_file_paths_targeted(session_path, 'postprobe', 'LED');
    postprobe_sessions = 1:2;
    kept_postprobes = [];
    for is = postprobe_sessions
        kept_postprobes = [kept_postprobes; all_postprobes(contains(all_postprobes, ['probe_0' num2str(is)]))];
    end

    kept_problems_first = [];
    kept_problems_last = [];
    problem_sessions = 1:6;
    for is = problem_sessions
        all_problems = [get_file_paths_targeted(session_path, ['mevar0' num2str(is)], 'mat')...
            get_file_paths_targeted(session_path, ['hivar0' num2str(is)], 'mat')];
        kept_problems_first = [kept_problems_first; all_problems(1)];
        kept_problems_last = [kept_problems_last; all_problems(end)];
    end
session_list = [kept_preprobes; kept_postprobes; kept_problems_first; kept_problems_last];
session_list = session_list([1 9 15 2 10 16 3 11 17 4 12 18 5 13 19 6 14 20 7 8])
cell_reg_mtx = cell_reg_mtx(:, [1 9 15 2 10 16 3 11 17 4 12 18 5 13 19 6 14 20 7 8]);
cell_reg_mtx = cell_reg_mtx(:, sessions);
session_list = session_list(sessions);


%}


%% Compute and activity matrices

% preallocate
session_mtx_cell = cell(size(session_list,1), 1);

% iterate through sessions
for isesh = 1:size(session_list,1)
    
    % trial activity
    trial_activity_mtx = image_trial_activity_mtx(session_list{isesh});
    
    % remove empty trials
    del_trl_idx = false(size(trial_activity_mtx,3),1);
    for itrl = 1:size(trial_activity_mtx,3)
        trl_vect = trial_activity_mtx(:,:,itrl); 
        if all(isnan(trl_vect(:)))
            del_trl_idx(itrl) = true;
        end
    end
    trial_activity_mtx(:,:,del_trl_idx) = [];
    
    % only first 41 trials
    trial_activity_mtx = trial_activity_mtx(:, :, 1:41);
    
    % remove time bins that are empty on any trial
    trial_activity_mtx = trial_activity_mtx(:, isfinite(sum(sum(trial_activity_mtx,3),1)), :);
    
    % only fixed period
    stage_bins = cumsum(time_series_event_spacing).*100;
    trial_activity_mtx = trial_activity_mtx(:, 1:stage_bins(4), :);
    
    % load
    session_mtx_cell{isesh} = trial_activity_mtx;
    
end



%% Reorganize sessions to match universal cell registry

% iterate through sessions
for isesh = 1:size(session_list,1)

    % session matrix with all cells
    session_mtx_cell_hold = zeros(size(cell_reg_mtx,1), size(trial_activity_mtx,2), size(trial_activity_mtx,3)); % matrix of zeros for inactive cells
    session_mtx_cell_hold(cell_reg_mtx(:,isesh)>0,:,:) = session_mtx_cell{isesh}(cell_reg_mtx(cell_reg_mtx(:,isesh)>0,isesh),:,:); % load active cells
    
    % load updated session matrix
    session_mtx_cell{isesh} = session_mtx_cell_hold;
    
end



%% Form indices for plotting

% preallocate session-level indices
session_number_idx = cell(size(session_mtx_cell));
trial_number_idx = cell(size(session_mtx_cell));
time_bin_idx = cell(size(session_mtx_cell));

% trial idx
trial_idx = nan(1,1,41); 
trial_idx(1,1,:) = 1:41;
trial_idx = repmat(trial_idx, size(session_mtx_cell_hold,1), size(session_mtx_cell_hold,2), 1);

% time idx
time_idx = repmat(1:size(session_mtx_cell_hold,2), size(session_mtx_cell_hold,1), 1, size(session_mtx_cell_hold,3));

% iterate through all sessions
for isesh = 1:length(session_mtx_cell)
    
    % preallocate cell-level indices
    session_number_idx{isesh} = repmat(isesh, size(session_mtx_cell{isesh}));
    trial_number_idx{isesh} = trial_idx;
    time_bin_idx{isesh} = time_idx;

end



%% Merge activity matrices

% merge neural activity
activity_mtx = cell_to_mtx(session_mtx_cell);

% merge session index
session_number_idx = cell_to_mtx(session_number_idx);
session_number_idx = session_number_idx(:,1);

% merge trial index
trial_number_idx = cell_to_mtx(trial_number_idx);
trial_number_idx = trial_number_idx(:,1);

% merge time index
time_bin_idx = cell_to_mtx(time_bin_idx);
time_bin_idx = time_bin_idx(:,1);



%% compute distances



end



%% internal functions

function mtx_out = cell_to_mtx(cell_in)
% function for merging the nested cells into a single matrix 
% (trials x concatenated neurons; see also indices)

mtx_out = [];
for isesh = 1:length(cell_in)
    for itrl = 1:size(cell_in{isesh},3)
        trl_hold = cell_in{isesh}(:,:,itrl);
        mtx_out = [mtx_out; trl_hold(:)'];
    end
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










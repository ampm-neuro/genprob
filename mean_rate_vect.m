function [mean_rate_matrix] = mean_rate_vect(mevar_or_hivar, subject_ids)
% computes the average firing rate vector of all cells in each stage

% for time warping
time_series_event_spacing = [0.2 1.1 1.0 2.0 8.5 2.0];
session_chron_reorder_idx = [1 9 15 2 10 16 3 11 17 4 12 18 5 13 19 6 14 20 7 8];



%% get all session paths
all_paths = cell(1,length(subject_ids));
for isubj = 1:length(subject_ids)
    
    % get session paths
    last_prob_sessions = [];
    for iprob = 1:6
        prob_sessions = get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_' mevar_or_hivar '_imaging_hpc\' subject_ids{isubj} ], ['var0' num2str(iprob)], 'LED');
        last_prob_sessions = [last_prob_sessions; prob_sessions(end)];
    end
    session_cell = [...
        get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_' mevar_or_hivar '_imaging_hpc\' subject_ids{isubj} ], 'preprobe', 'LED');...
        get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_' mevar_or_hivar '_imaging_hpc\' subject_ids{isubj} ] , 'postprobe', 'LED');...
        get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_' mevar_or_hivar '_imaging_hpc\' subject_ids{isubj} ], 'var0', '-01', 'LED');...
        last_prob_sessions];
    
    % load chronologically
    if length(session_cell)==19
        scri = session_chron_reorder_idx(1:end-1); scri(scri>=8) = scri(scri>=8)-1;
        all_paths{isubj} = session_cell(scri);
    else
        all_paths{isubj} = session_cell(session_chron_reorder_idx);
    end
    
    all_paths{isubj}
    
end

%% iterate through sessions compute average rate vect
mean_rate_matrix = [];
for istage = 1:20
    stage_rate_vects = [];
    for isubj = 1:length(subject_ids)
        
        % check if subj finished session
        if length(all_paths{isubj})<istage
            continue
        end
        
        % load session
        load(all_paths{isubj}{istage})

        % session matrix
        [sorted_session_mtx] = image_mean_activity_timewarp(...
        trl_mtx, trl_idx, frame_times, traces, 1:size(traces,1), unique(trl_idx), time_series_event_spacing);
        
        stage_rate_vects = [stage_rate_vects; sorted_session_mtx];
    end
    mean_rate_matrix = [mean_rate_matrix; nanmean(stage_rate_vects)];
end
    
    
%% plots

%figure;
%ampm_pcolor(mean_rate_matrix)





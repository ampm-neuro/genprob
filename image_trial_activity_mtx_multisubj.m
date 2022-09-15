function [combine_trial_activity_mtx, subject_ids, session_nums, cell_regist_mtx_cell] = image_trial_activity_mtx_multisubj(subjects, sessions, trials, tses)
% analogous to image_trial_activity_mtx, but with concatenated sessions

event_frame = cumsum(tses(1:end-1)).*100;
session_chron_reorder_idx = [1 9 15 2 10 16 3 11 17 4 12 18 5 13 19 6 14 20 7 8];

% iterate through subjects
combine_trial_activity_mtx = [];
subject_ids = [];
session_nums = [];
cell_regist_mtx_cell = cell(1, length(subjects));
for isubj = 1:length(subjects)
    isubj
    
    % get session paths of interest
    %
        % load cell registration matrix
        load(['cellreg_data\' subjects{isubj} '\cell_reg_' subjects{isubj} '.mat'])
        cell_regist_mtx = cell_registered_struct.cell_to_index_map;

        % chron
        if strcmp(subjects{isubj}, '690330m1')
            cell_regist_mtx = [cell_regist_mtx(:,1:9) nan(size(cell_regist_mtx(:,1))) cell_regist_mtx(:,10:end)]; 
            scri = session_chron_reorder_idx;
        elseif size(cell_regist_mtx,2)==19
            scri = session_chron_reorder_idx(1:end-1); scri(scri>=8) = scri(scri>=8)-1;
        else
            scri = session_chron_reorder_idx;
        end

        % desired sessions only
        scri = scri(sessions);
        cell_regist_mtx_cell{isubj} = cell_regist_mtx(:, scri);

        % get session paths
        last_prob_sessions = [];
        for iprob = 1:6
            prob_sessions = [get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_imaging_hpc\' subjects{isubj} ], ['var0' num2str(iprob)], 'LED');...
                             get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_hivar_imaging_hpc\' subjects{isubj} ], ['var0' num2str(iprob)], 'LED');];
            last_prob_sessions = [last_prob_sessions; prob_sessions(end)];
        end
        session_cell = [...
            get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_imaging_hpc\' subjects{isubj} ], 'preprobe', 'LED');...
            get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_hivar_imaging_hpc\' subjects{isubj} ], 'preprobe', 'LED');...
            get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_imaging_hpc\' subjects{isubj} ] , 'postprobe', 'LED');...
            get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_hivar_imaging_hpc\' subjects{isubj} ] , 'postprobe', 'LED');...
            get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_imaging_hpc\' subjects{isubj} ], 'var0', '-01', 'LED');...
            get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_hivar_imaging_hpc\' subjects{isubj} ], 'var0', '-01', 'LED');...
            last_prob_sessions];
        session_cell = session_cell(scri);
        
        
    % iterate through sessions
    %
        for isesh = 1:size(session_cell,1)
            
            % compute activity
            trial_activity_mtx = image_trial_activity_mtx(session_cell{isesh}, tses);
            
                % only trials and time bins (np2 to end of random delay) of interest 
                trial_activity_mtx = trial_activity_mtx(:, sum(~isnan(nanmean(trial_activity_mtx,3)))>1, trials);
                
                trial_activity_mtx = trial_activity_mtx(:, event_frame(1):event_frame(end-1), :);
            
            % load
            combine_trial_activity_mtx = [combine_trial_activity_mtx; trial_activity_mtx];
            
            % update indices
            subject_ids = [subject_ids; repmat(isubj, size(trial_activity_mtx,1), 1)];
            session_nums = [session_nums; repmat(sessions(isesh), size(trial_activity_mtx,1), 1)];
        
        end
end     

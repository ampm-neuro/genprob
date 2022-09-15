function [within_session_corr_means, cell_session_means] = within_session_reliability_frfields(subject_ids, sessions, tses)
% compute pairwise correlations between every trial within a session
% correlations are the average of each cells individual correlation of its
% activity on each trial with its activity on every other trial

event_frame = cumsum([tses(1:end-1)]).*100;
session_chron_reorder_idx = [1 9 15 2 10 16 3 11 17 4 12 18 5 13 19 6 14 20 7 8];

cell_session_means = cell(length(subject_ids), length(sessions));
within_session_corr_means = nan(length(subject_ids), length(sessions));
for isubj = 1:length(subject_ids)
    subject_ids{isubj}
    
    % get session list
    %
    
        % load cell registration matrix
        load(['cellreg_data\' subject_ids{isubj} '\cell_reg_' subject_ids{isubj} '.mat'])
        cell_regist_mtx = cell_registered_struct.cell_to_index_map;

        % chron
        if strcmp(subject_ids{isubj}, '690330m1')
            cell_regist_mtx = [cell_regist_mtx(:,1:9) nan(size(cell_regist_mtx(:,1))) cell_regist_mtx(:,10:end)]; 
            scri = session_chron_reorder_idx;
        elseif size(cell_regist_mtx,2)==19
            scri = session_chron_reorder_idx(1:end-1); scri(scri>=8) = scri(scri>=8)-1;
        else
            scri = session_chron_reorder_idx;
        end

        % desired sessions only
        scri = scri(sessions);
        cell_regist_mtx = cell_regist_mtx(:, scri);

        % get session paths
        last_prob_sessions = [];
        for iprob = 1:6
            prob_sessions = [get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_imaging_hpc\' subject_ids{isubj} ], ['var0' num2str(iprob)], 'LED');...
                             get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_hivar_imaging_hpc\' subject_ids{isubj} ], ['var0' num2str(iprob)], 'LED');];
            last_prob_sessions = [last_prob_sessions; prob_sessions(end)];
        end
        session_cell = [...
            get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_imaging_hpc\' subject_ids{isubj} ], 'preprobe', 'LED');...
            get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_hivar_imaging_hpc\' subject_ids{isubj} ], 'preprobe', 'LED');...
            get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_imaging_hpc\' subject_ids{isubj} ] , 'postprobe', 'LED');...
            get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_hivar_imaging_hpc\' subject_ids{isubj} ] , 'postprobe', 'LED');...
            get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_imaging_hpc\' subject_ids{isubj} ], 'var0', '-01', 'LED');...
            get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_hivar_imaging_hpc\' subject_ids{isubj} ], 'var0', '-01', 'LED');...
            last_prob_sessions];
        session_cell = session_cell(scri);
        
     
     % compute mean pairwise correlation for each session
     %
     for isesh = 1:size(session_cell,1)
                  
         % load
         load(session_cell{isesh})
         
         % compute activity of every cell on every trial
         hm_cell = tw_activity_plot_trial_hm_warp(trl_mtx, trl_idx, frame_times, traces, 1:size(traces,1), 1:size(trl_mtx,1), ceil([sum(tses(1:3)) sum(tses(4:end))]), tses);
         
         % compute cell means
         cell_trial_means = nan(((size(hm_cell{1},1)-1)^2 + (size(hm_cell{1},1)-1))/2, length(hm_cell));
         for ic = 1:length(hm_cell)
            
             % activity during period from second np to end of random delay
             hm = hm_cell{ic}(:, sum(~isnan(hm_cell{ic}))>1);
             hm = hm(:, event_frame(1):event_frame(end-1));
             
             % pairwise over all trials
             trial_comp = 0;
             for itrials1 = 1:size(hm_cell{ic},1)
                 for itrials2 = 1:size(hm_cell{ic},1)
                    % unique only
                    if itrials2<=itrials1
                         continue
                    else
                        trial_comp = trial_comp + 1;
                    end
                    
                    % load trial mean for this cell

                    nnan_idx = ~isnan(hm(itrials1,:)) & ~isnan(hm(itrials2,:));
                    if sum(nnan_idx)==0
                        continue
                    end
                    cell_trial_means(trial_comp, ic) = corr(hm(itrials1,nnan_idx)', hm(itrials2,nnan_idx)');
                     
                 end
             end

         end
         cell_session_means{isubj, isesh} = nanmean(cell_trial_means)';
         within_session_corr_means(isubj, isesh) = nanmean(nanmean(cell_trial_means));
     end
    
    
    
    
    
end




end


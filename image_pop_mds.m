function [distance_mtx, session_mtx_cell, session_id, trial_id] = image_pop_mds(subject_id, cell_reg_mtx, sessions, session_list, time_series_event_spacing)
% plot distance matrix and MDS for every trial from a single subject



%% end of fixed delay
eod = cumsum(time_series_event_spacing);
eod = eod(4)/eod(end);


%% get sessions for subject
%
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
all_problems_first = get_file_paths_targeted(session_path, 'mevar', '01.mat');
    problem_sessions = sessions(sessions>=9)-8;
    kept_problems_first = [];
    for is = problem_sessions
        kept_problems_first = [kept_problems_first; all_problems_first(contains(all_problems_first, ['mevar0' num2str(is)]))];
    end
kept_problems_last = [];
for iprob = 1:6
prob_sessions = get_file_paths_targeted(session_path, ['var0' num2str(iprob)], 'LED');
kept_problems_last = [kept_problems_last; prob_sessions(end)];
end

kept_problems_all = sort([kept_problems_first; kept_problems_last]);
    
session_list = [kept_preprobes; kept_postprobes; kept_problems_all]
%}


%% Compute activity matrices
%
% preallocate
session_mtx_cell = cell(size(session_list,1), 1);

% iterate through sessions
for isesh = 1:size(session_list,1)
    
    % load session
    load(session_list{isesh})

    % only reward-unavailable trials?
    hm_cell_trials = unique(trl_idx);
    %hm_cell_trials = intersect(hm_cell_trials, find(trl_mtx(:,3)==0));
    
    % session matrix for every cell
    hm_cells = tw_activity_trial_hm_full_warp(trl_mtx, trl_idx, frame_times, traces, 1:size(traces,1), hm_cell_trials, [3 12], time_series_event_spacing);
    
    % remove times bins with nans on any trial
    for icell = 1:length(hm_cells)
        %hm_cells{icell} = hm_cells{icell}(:,~isnan(sum(hm_cells{icell},1)));
        %hm_cells{icell} = reshape(zscore_mtx(hm_cells{icell}(:)), size(hm_cells{icell}));
    end

    % load
    session_mtx_cell{isesh} = hm_cells;
    
end
%save(['image_pop_mds_full_save_' subject_id '.mat'], 'session_mtx_cell')
%error
%}
%load(['image_pop_mds_full_save_' subject_id '.mat'], 'session_mtx_cell')

%% Reorganize sessions to match universal cell registry
%
% iterate through sessions
for isesh = 1:size(session_list,1)
    
    % matrix of zeros for inactive cells
    inactive_mtx = zeros(size(session_mtx_cell{isesh}{1}));

    % session matrix with all cells
    session_mtx_cell_hold = cell(1, size(cell_reg_mtx,1));
    
    % preload all cells as inactive
    for icell = 1:length(cell_reg_mtx(:,sessions(isesh)))
        session_mtx_cell_hold{icell} = inactive_mtx;
    end
    
    % load active cells
    
    isesh
    
    for icell = 1:size(size(session_mtx_cell{isesh}),2) %unique(cell_reg_mtx(cell_reg_mtx(:,sessions(isesh))>0,sessions(isesh)))'
        
        icell
        size(session_mtx_cell)
        size(session_mtx_cell{isesh})
    	session_mtx_cell_hold{find(cell_reg_mtx(:,sessions(isesh))==icell, 1)} = session_mtx_cell{isesh}{icell};
    end
    
    % load updated session matrix
    session_mtx_cell{isesh} = session_mtx_cell_hold;
    
end
['image_pop_mds_full_save_' subject_id '.mat']
save(['image_pop_mds_full_save_' subject_id '.mat'], 'session_mtx_cell', '-v7.3')
error
%}
%load(['image_pop_mds_full_save_' subject_id '.mat'], 'session_mtx_cell')

%{
for isesh = 1:length(session_mtx_cell)
    size(session_mtx_cell{isesh}{1})
end
%}


%% Average activity 

% session and trial counts
session_id = [];
trial_id = [];

min_trials = 41;

% iterate through sessions
for isesh = 1:length(session_mtx_cell)
    eod_sesh = ceil(eod*size(session_mtx_cell{isesh}{1},2));
    
    % iterate through cells
    for icell = 1:length(session_mtx_cell{isesh})
        
        % all activity during fixed delay, first 41 trials
        session_activity_earlyonly = session_mtx_cell{isesh}{icell}(:,1:eod_sesh);
        session_activity_first41 = session_activity_earlyonly(1:min_trials,:);
        session_mtx_cell{isesh}{icell} = zscore_mtx(session_activity_first41')';

    end

    % each session is a 3d matrix of activity: (trials, time, cells)
    session_mtx_cell{isesh} = ...
        reshape(cell2mat(session_mtx_cell{isesh}), size(session_mtx_cell{isesh}{1},1), size(session_mtx_cell{isesh}{1},2), length(session_mtx_cell{isesh}));
   
    % indices
    session_id = [session_id; repmat(isesh, size(session_mtx_cell{isesh},1), 1)];
    trial_id = [trial_id; (1:size(session_mtx_cell{isesh},1))'];
end



%% Compute distance between every trial in every session
%
% preallocate
num_trials = size(trial_id,1);
distance_mtx = zeros(size(trial_id,1));

% iterate through all session trials
for isesh1 = 1:size(session_list,1)
    session1 = session_mtx_cell{isesh1};
    
    disp(num2str(isesh1))
    
    for isesh2 = 1:size(session_list,1)

            if isesh1 < isesh2
                continue
            end

            session2 = session_mtx_cell{isesh2};
    
        trl1_ct = sum(session_id<isesh1);
        for itrl1 = 1:size(session1,1)
            trl1_ct = trl1_ct+1;

            % (cells, time)
            trial1 = squeeze(session1(itrl1,:,:))';

            trl2_ct = sum(session_id<isesh2);
            for itrl2 = 1:size(session2,1)
                trl2_ct = trl2_ct+1;
                
                if isesh1 == isesh2 && itrl1 <= itrl2
                    continue
                end
                
                % (cells, time)
                trial2 = squeeze(session2(itrl2,:,:))';
                
                % common cells only
                common_cell_idx = ~isnan(trial1(:,1)) & ~isnan(trial2(:,1));
                trial1_common = trial1(common_cell_idx,:);
                trial2_common = trial2(common_cell_idx,:);
                
                % compute difference
                trial_diff_common = trial1_common-trial2_common;
                trial_diff_common  = trial_diff_common.^2;
                trial_diff_common = mean(trial_diff_common(:));
                
                % load difference                
                %distance_mtx(trl1_ct, trl2_ct) = sqrt(sum(sum((trial1_common-trial2_common).^2)))./sqrt(numel(trial1_common));
                distance_mtx(trl1_ct, trl2_ct) = trial_diff_common; %./sqrt(numel(trial1_common));
                
                %ampm_pcolor(distance_mtx); drawnow
                
            end
        end
    end
end

% fill in givens
distance_mtx = sum(cat(3, distance_mtx, distance_mtx'),3);

save(['image_pop_mds_dm_save_' subject_id '.mat'], 'distance_mtx', 'session_id')
%}
%load(['image_pop_mds_dm_save_' subject_id '.mat'], 'distance_mtx', 'session_id')


%% plot matrix

figure; 
ampm_pcolor(distance_mtx); 
axis square
caxis([0 2.5])
colorbar
hold on
for xpos = 1:41:41*21
    plot(xpos.*[1 1], ylim, 'k-')
end
for ypos = 1:41:41*21
    plot(xlim, ypos.*[1 1], 'k-')
end
drawnow



%% plot mds

[mds_coords, stress] = mds_plot_2(distance_mtx, 3, session_id);











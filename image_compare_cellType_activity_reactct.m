function [dispersion_mtx, reliability_mtx, reactivation_ct_mtx] = image_compare_cellType_activity_reactct(session_list, cell_regist_mtx, time_series_event_spacing)
% plot the information properties of topDown and bottomUp cells
% consider sorting inputs chronologically



%% compute trial activity for every cell in every session
%
% preallocate
session_mtx_cell = cell(size(session_list,1), 1);

% iterate through sessions
for isesh = 1:size(session_list,1)
    
    % load session
    load(session_list{isesh})

    % all trials
    hm_cell_trials = unique(trl_idx);
    
    % session matrix for every cell
    hm_cells = tw_activity_trial_hm_full_warp(trl_mtx, trl_idx, frame_times, traces, 1:size(traces,1), hm_cell_trials, [3 12], time_series_event_spacing);
    
    % remove times bins with nans on any trial
    for icell = 1:length(hm_cells)
        hm_cells{icell} = hm_cells{icell}(:,~isnan(sum(hm_cells{icell},1)));
    end

    % load
    session_mtx_cell{isesh} = hm_cells;
    
end
%save('comp_cell_type_info.mat', '-v7.3')
disp('1')
%}



%% Reorganize sessions to match universal cell registry
%
% iterate through sessions
sessions = 1:size(session_list,1);
for isesh = 1:size(session_list,1)
    
    % matrix of zeros for inactive cells
    inactive_mtx = zeros(size(session_mtx_cell{isesh}{1}));

    % session matrix with all cells
    session_mtx_cell_hold = cell(size(cell_regist_mtx,1),1);
    
    % preload all cells as inactive
    for icell = 1:length(cell_regist_mtx(:,sessions(isesh)))
        session_mtx_cell_hold{icell} = inactive_mtx;
    end
    
    % load active cells
    for icell = unique(cell_regist_mtx(cell_regist_mtx(:,sessions(isesh))>0,sessions(isesh)))'
    	session_mtx_cell_hold{find(cell_regist_mtx(:,sessions(isesh))==icell, 1)} = session_mtx_cell{isesh}{icell};
    end
    
    % load updated session matrix
    session_mtx_cell{isesh} = session_mtx_cell_hold;
    
end

%save('comp_cell_type_info.mat', '-v7.3')
disp('1.5')
%}



%% compute dispersion and reliability
%
% dispersion, not info content
[dispersion_mtx, reliability_mtx] = all_cells_infoContent(1:size(cell_regist_mtx,1), session_mtx_cell);

% set dispersion to 0 : 1
num_time_bins = 1279; %size(session_mtx_cell{1}{1},2);
divisor = 10;
% (1-proportion of the contiguous time bins required to hold 50% of activity)
dispersion_mtx = 1-dispersion_mtx./ceil(num_time_bins/divisor);

%save('comp_cell_type_info.mat', '-v7.3')
disp('2')
%}


%% label cells by number of reactivations
reactivation_ct_mtx = cell_regist_mtx;
for ic = 1:size(reactivation_ct_mtx,1)
    reactivation_ct_mtx(ic, reactivation_ct_mtx(ic,:)>0) = 1:sum(reactivation_ct_mtx(ic,:)>0);
end







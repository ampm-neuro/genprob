function [common_cell_corrs, common_cell_corrs_means] = cell_turnover_timewarp_trials_cellTrajectories(session_cell, cell_registration, time_series_event_spacing)
% plots average timewarped trials for every registered cell in every
% session, with inactive cells receiving zeros for activity

session_cell
probe_nums = [1:3:19 20];
problem_nums = setdiff(1:20, probe_nums);

%% compute session matrix for each session
%{
session_mtx_cell = cell(size(session_cell));
session_max_cell = cell(size(session_cell));
session_min_cell = cell(size(session_cell));
mean_activity_mtx = [];
for isesh = 1:length(session_cell)
   
    % load session
    load(session_cell{isesh})
    
    % session matrix
    [sorted_session_mtx, srt_idx, ~, max_time_bin_unsorted, ~, min_time_bin_unsorted] = image_mean_activity_timewarp(...
    trl_mtx, trl_idx, frame_times, traces, 1:size(traces,1), unique(trl_idx), time_series_event_spacing, isesh);

    % unsort
    unsort_idx = nan(size(srt_idx));
    unsort_idx(srt_idx) = 1:length(srt_idx);
    session_mtx_cell{isesh} = sorted_session_mtx(unsort_idx,:);
    session_max_cell{isesh} = max_time_bin_unsorted;
    session_min_cell{isesh} = min_time_bin_unsorted;
    
    mean_activity_mtx=[mean_activity_mtx ; mean(sorted_session_mtx)];
    
end

% average cell activity
figure; hold on
ampm_pcolor(mean_activity_mtx); title('0')
% red line
event_frame = cumsum([0 time_series_event_spacing]).*100;
for i = event_frame 
   plot(i.*[1 1], ylim, 'r-', 'linewidth', 1.5) 
end
xlim([0.5 size(all_mtx,2)+0.5])
%error


save('cellTrajectories.mat', 'session_mtx_cell', 'session_max_cell', 'session_min_cell')
%}
load('cellTrajectories.mat', 'session_mtx_cell', 'session_max_cell', 'session_min_cell')



%% plot each cell's firing rate maps
%

% iterate through cells
[~,cellnums_sorted] = sort(sum(cell_registration>0,2), 'descend');
%persist_to_punctate = [274 237 243 168 15 230 64 228 141 224];
for icell = 1:size(cell_registration,1) %persist_to_punctate

    figure; hold on
    cell_history_mtx = nan(length(session_cell), size(session_mtx_cell{1},2));    
    
    %figure; hold on
    % iterate through sessions 
    for isesh = 1:20
    
        if cell_registration(icell,isesh) == 0
            continue
        end
        
        cell_history_mtx(isesh, :) = session_mtx_cell{isesh}(cell_registration(icell,isesh), :);

    end
    
    % plot
    imagesc(cell_history_mtx);

    % red line
    event_frame = cumsum([0 time_series_event_spacing]).*100;
    xticks_hold = [];
    for i = event_frame 
       plot(i.*[1 1], ylim, 'r-', 'linewidth', 1.0) 
       xticks_hold = [xticks_hold i];
    end

    % aesthetics
    set(gca, 'YDir','reverse')
    ylim([0+0.5 size(cell_history_mtx,1)+0.5])
    xlim([0.5 size(cell_history_mtx,2)+0.5])
    xticks(xticks_hold)
    xticklabels_universe = {'NP1', 'NP2', 'HE1', 'TO', 'RD', 'RWD', 'TE'};
    xticklabels(xticklabels_universe)
    yticks(probe_nums)
    yticklabels({'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8'})
    title(num2str(icell))

end


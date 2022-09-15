function [common_cell_corrs, common_cell_corrs_means] = cell_turnover_timewarp_trials_cellTrajectories_rp(session_cell, cell_registration, time_series_event_spacing)
% plots average timewarped trials for every registered cell in every
% session, with inactive cells receiving zeros for activity

session_cell

probe_nums = [1:3:19 20];
problem_nums = setdiff(1:20, probe_nums);

%% compute rich and poor session matrix cells
%{
session_mtx_cell_rich = cell(size(session_cell));
session_max_cell_rich = cell(size(session_cell));
session_min_cell_rich = cell(size(session_cell));
session_mtx_cell_poor = cell(size(session_cell));
session_max_cell_poor = cell(size(session_cell));
session_min_cell_poor = cell(size(session_cell));
for isesh = 1:length(session_cell)
   
    % load session
    load(session_cell{isesh})
    
    % rich and poor tones
    if contains(session_cell{isesh}, 'probe')
        rich_trials = unique(trl_idx);
        poor_trials = unique(trl_idx);
    else
        
        unq_tones = unique(floor(trl_mtx(:,2)));
        reinforced_trial_tone_counts = histcounts(trl_mtx(trl_mtx(:,3)==1,2), [unq_tones(1) unq_tones(2)+1]);
        rich_tone = unq_tones(reinforced_trial_tone_counts==max(reinforced_trial_tone_counts));
        poor_tone = setdiff(unq_tones, rich_tone);
        rich_trials = intersect(unique(trl_idx), find(floor(trl_mtx(:,2))==rich_tone));
        poor_trials = intersect(unique(trl_idx), find(floor(trl_mtx(:,2))==poor_tone));
    end
    
    
    % session matrix rich
    %
    
        [sorted_session_mtx, srt_idx, ~, max_time_bin_unsorted, ~, min_time_bin_unsorted] = image_mean_activity_timewarp(...
        trl_mtx, trl_idx, frame_times, traces, 1:size(traces,1), rich_trials, time_series_event_spacing);

        % unsort
        unsort_idx = nan(size(srt_idx));
        unsort_idx(srt_idx) = 1:length(srt_idx);
        session_mtx_cell_rich{isesh} = sorted_session_mtx(unsort_idx,:);
        session_max_cell_rich{isesh} = max_time_bin_unsorted;
        session_min_cell_rich{isesh} = min_time_bin_unsorted;
        
        
    % session matrix poor
    %
    
        [sorted_session_mtx, srt_idx, ~, max_time_bin_unsorted, ~, min_time_bin_unsorted] = image_mean_activity_timewarp(...
        trl_mtx, trl_idx, frame_times, traces, 1:size(traces,1), poor_trials, time_series_event_spacing);

        % unsort
        unsort_idx = nan(size(srt_idx));
        unsort_idx(srt_idx) = 1:length(srt_idx);
        session_mtx_cell_poor{isesh} = sorted_session_mtx(unsort_idx,:);
        session_max_cell_poor{isesh} = max_time_bin_unsorted;
        session_min_cell_poor{isesh} = min_time_bin_unsorted;
    
end
save('cellTrajectories_rp_nonnorm.mat', 'session_mtx_cell_rich', 'session_max_cell_rich', 'session_min_cell_rich', 'session_mtx_cell_poor', 'session_max_cell_poor', 'session_min_cell_poor')
%}
load('cellTrajectories_rp_nonnorm.mat', 'session_mtx_cell_rich', 'session_max_cell_rich', 'session_min_cell_rich', 'session_mtx_cell_poor', 'session_max_cell_poor', 'session_min_cell_poor')


%% plot each cell's firing rate maps
%

% iterate through cells

cell_nums_prob = find(sum(cell_registration(:, problem_nums),2)>0);
[~,cellnums_sorted] = sort(sum(cell_registration(cell_nums_prob,:)>0,2), 'descend');
cell_nums_prob = cell_nums_prob(cellnums_sorted);
persist_to_punctate = [274 237 243 168 15 230 64 228 141 224];
for icell =  cell_nums_prob(1:5:end)' % 1:size(cell_registration,1) %[6 8 9 32 35 83 96] %

    figure

    cell_history_mtx_rich = nan(length(session_cell), size(session_mtx_cell_rich{1},2));
    cell_history_mtx_poor = nan(length(session_cell), size(session_mtx_cell_poor{1},2));
    
    % iterate through sessions 
    for isesh = 1:20
    
        if cell_registration(icell,isesh) == 0
            continue
        elseif ismember(isesh, [1:3:19 20])
            continue
        end
        
        rich_vect = session_mtx_cell_rich{isesh}(cell_registration(icell,isesh), :);
        poor_vect = session_mtx_cell_poor{isesh}(cell_registration(icell,isesh), :);
        
        % normalize
        rich_poor_norm = [rich_vect poor_vect];
        rich_poor_norm = rich_poor_norm-min(rich_poor_norm);
        rich_poor_norm = rich_poor_norm./max(rich_poor_norm);
        cell_history_mtx_rich(isesh, :) = rich_poor_norm(1:length(rich_vect));
        cell_history_mtx_poor(isesh, :) = rich_poor_norm(length(rich_vect)+1 : end);
    end
    
    % plot
    % rich trials
    %
    subplot(1,2,1); hold on; title rich
    imagesc(cell_history_mtx_rich);

    % red line
    event_frame = cumsum([0 time_series_event_spacing]).*100;
    xticks_hold = [];
    for i = event_frame 
       plot(i.*[1 1], ylim, 'r-', 'linewidth', 1.0) 
       xticks_hold = [xticks_hold i];
    end

    % aesthetics
    set(gca, 'YDir','reverse')
    ylim([0+0.5 size(cell_history_mtx_rich,1)+0.5])
    xlim([0.5 size(cell_history_mtx_rich,2)+0.5])
    xticks(xticks_hold)
    xticklabels_universe = {'NP1', 'NP2', 'HE1', 'TO', 'RD', 'RWD', 'TE'};
    xticklabels(xticklabels_universe)    
    caxis_hold_rich = caxis;
    
    
    
    % plot
    % poor trials
    %
    subplot(1,2,2); hold on; title poor
    imagesc(cell_history_mtx_poor);

    % red line
    event_frame = cumsum([0 time_series_event_spacing]).*100;
    xticks_hold = [];
    for i = event_frame 
       plot(i.*[1 1], ylim, 'r-', 'linewidth', 1.0) 
       xticks_hold = [xticks_hold i];
    end

    % aesthetics
    set(gca, 'YDir','reverse')
    ylim([0+0.5 size(cell_history_mtx_poor,1)+0.5])
    xlim([0.5 size(cell_history_mtx_poor,2)+0.5])
    xticks(xticks_hold)
    xticklabels_universe = {'NP1', 'NP2', 'HE1', 'TO', 'RD', 'RWD', 'TE'};
    xticklabels(xticklabels_universe)
    yticks(probe_nums)
    yticklabels({'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8'})
    caxis_hold_poor = caxis;
    
    % meta title
    sgtitle(num2str(icell))
    set(gcf, 'Position', [528 766 1032 572])
    
    % common caxis
    caxis_hold = [caxis_hold_rich(:); caxis_hold_poor(:)];
    caxis_hold = [min(caxis_hold) max(caxis_hold)];
    subplot(1,2,1); caxis(caxis_hold)
    subplot(1,2,2); caxis(caxis_hold)

end


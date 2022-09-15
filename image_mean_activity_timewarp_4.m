function image_mean_activity_timewarp_4(trl_mtx, trl_idx, frame_times, traces, neurons)
% runs image_mean_activity_timewarp for each trial type

% indices

%intersect(unique(trl_idx), find(rich_trl_idx(trl_mtx)) find(trl_mtx(:,3)==1))

rich_idx =  rich_trl_idx(trl_mtx);
rwd_idx = ~isnan(trl_mtx(:,11));
probe_idx = ~isnan(trl_mtx(:,11));

% rich, reward
trials = find(rich_idx & rwd_idx); trials = trials(ismember(trials, unique(trl_idx)));
[~, sort_idx, peaks] = image_mean_activity_timewarp(trl_mtx, trl_idx, frame_times, traces,...
    neurons, trials);
title('rich, reward')

%{
% rich no reward
trials = find(rich_idx & ~rwd_idx); trials = trials(ismember(trials, unique(trl_idx)));
image_mean_activity_timewarp(trl_mtx, trl_idx, frame_times, traces,...
    neurons, trials, sort_idx);
title('rich, no reward')
%


% poor, reward
trials = find(~rich_idx & rwd_idx); trials = trials(ismember(trials, unique(trl_idx)));
image_mean_activity_timewarp(trl_mtx, trl_idx, frame_times, traces,...
    neurons, trials, sort_idx);
title('poor, reward')

%
% poor, no reward
trials = find(~rich_idx & ~rwd_idx); trials = trials(ismember(trials, unique(trl_idx)));
image_mean_activity_timewarp(trl_mtx, trl_idx, frame_times, traces,...
    neurons, trials, sort_idx);
title('poor, no reward')
%}


% correlations
tw_activity_corr_waits_warp_richpoor(trl_mtx, trl_idx, frame_times, traces, 1:size(traces,1), 3, [5 12], [0.5 1.1 1.5 2.0 3.5 2.00], sort_idx, peaks);
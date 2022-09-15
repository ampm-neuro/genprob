function image_mean_activity_timewarp_byprobe(trl_mtx, trl_idx, frame_times, traces, neurons)
% runs image_mean_activity_timewarp for trials when the animal waited
% above or below average, with or without rwd

% wait times
[wait_times, frequencies] = wait_times_prep(trl_mtx, 1);
mwt = mean(wait_times);

% indices
long_wait_freqs = frequencies(wait_times>mwt)
longpwait_idx =  ismember(floor(trl_mtx(:,2)), long_wait_freqs);
rwd_idx = ~isnan(trl_mtx(:,11));

% long probe wait, reward
[~, sort_idx] = image_mean_activity_timewarp(trl_mtx, trl_idx, frame_times, traces,...
    neurons, find(longpwait_idx & rwd_idx));
title('long-wait, reward')

% long probe wait, no reward
image_mean_activity_timewarp(trl_mtx, trl_idx, frame_times, traces,...
    neurons, find(longpwait_idx & ~rwd_idx), sort_idx);
title('long-wait, no reward')

%{
% short probe wait, reward
image_mean_activity_timewarp(trl_mtx, trl_idx, frame_times, traces,...
    neurons, find(~longpwait_idx & rwd_idx), sort_idx);
title('short-wait, reward')

% short probe wait, no reward
image_mean_activity_timewarp(trl_mtx, trl_idx, frame_times, traces,...
    neurons, find(~longpwait_idx & ~rwd_idx), sort_idx);
title('short-wait, no reward')
%}

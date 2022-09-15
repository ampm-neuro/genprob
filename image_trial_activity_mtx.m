function trial_activity_mtx = image_trial_activity_mtx(session_path, tses)
% outputs cell X timebin X trials

% load from session file
load(session_path, 'traces', 'frame_times', 'trl_idx', 'trl_mtx')

%timewarp traces
[traces, frame_times, trl_idx, trl_mtx] = timewarp_traces(traces, frame_times, trl_idx, trl_mtx, tses);

%activity matrix
neurons = 1:size(traces,1);
trials = unique(trl_idx);
tw_bounds =  [sum(tses(1:2))+5 sum(tses(3:end))+5];
%act_mtx = tw_activity_II(trl_mtx, trl_idx, frame_times, traces, neurons, trials, 3, tw_bounds);
trial_activity_cell = tw_activity_trial_hm_full_warp(trl_mtx, trl_idx, frame_times, traces, neurons, trials, tw_bounds, tses);

% cell 2 mat
trial_activity_mtx = nan(size(trial_activity_cell{1},1), size(trial_activity_cell{1},2), length(trial_activity_cell));
for itrial = 1:length(trial_activity_cell)
    trial_activity_mtx(:,:,itrial) = trial_activity_cell{itrial};
end

% trials, time, neuron
trial_activity_mtx = permute(trial_activity_mtx,[3,2,1]);
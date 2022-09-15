function hm_cell = tw_activity_trial_hm_warp(trl_mtx, trl_idx, all_pkl_frames, traces, neurons, trials, tw_bounds)
% runs tw_activity and produces one plot for each neuron comparing activity
% for each unique tone
%
% align events:
% NP onset = 1
% NP offset = 2
% Tone on = 3
% Head entry = 4
% Reward delivery = 5
% Reward receipt = 6
% Trial end = 7
%
% tw_bounds is a two item vector. the first item is the number of seconds
% before the align event the tw should start and the second item is the
% number of seconds after the align event the tw should end. eg [2 3] means
% the tw starts 2 seconds and ends 3s after the align event.
% compute all activity
%
% input tses is time series event spacing used to warp traces
%

%activity matrix
act_mtx = tw_activity_II(trl_mtx, trl_idx, all_pkl_frames, traces, neurons, trials, 4, tw_bounds);

% iterate through each neuron
for ic = 1:length(neurons)
    
    % iterate through tones
    hm_cell = nan(length(trials),size(act_mtx,2));
    for itrl = 1:length(trials)
    
        % compute trace
        trace_local = act_mtx(ic, :, itrl);
        % trace_local = norm_mtx(trace_local);
        hm_cell(itrl,:) = trace_local;
    end
end

    
function act_mtx = tw_activity(trl_mtx, trl_idx, medass_cell, all_pkl_frames, traces, neurons, trials, align_event, tw_bounds)
% computes a matrix of activities over time defined by tw_bounds
%
% INPUT:
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
%
% output:
% act_mtx = (neuron, activity over time, trial)


% remove abberations
%{
for i = 1:size(traces,1)
    traces(i,:) = remove_abberations(traces(i,:), 2.5);
end
%}

% all cells
all_cells = neurons;

% normalize all traces to [0 1]
%traces = traces(neurons,:);
traces = traces - min(traces, [], 2);
traces = traces ./ max(traces, [], 2);

% preallocate output matrix (activity sample every 0.01s)
num_samps = sum(tw_bounds)/0.01;
act_mtx = nan(length(neurons), num_samps, length(trials));

% iterate through all neurons
for ic = 1:length(neurons)
    
    current_cell = all_cells(ic);
    
    % iterate through all trials
    for itrl = 1:length(trials)
    
        % current items
        current_trial = trials(itrl);
        
        % time of align event
        switch align_event
            case 1 % NP onset
                ae_time = sum(trl_mtx(current_trial, [1 6]));
            case 2 % NP offset
                ae_time = sum(trl_mtx(current_trial, [1 9]));
            case 3 % Tone on
                ae_time = sum(trl_mtx(current_trial, [1 7]));
            case 4 % Head entry
                ae_time = sum(trl_mtx(current_trial, [1 10]));
            case 5 % Reward delivery
                ae_time = sum(trl_mtx(current_trial, [1 11]));
            case 6 % Reward receipt (first lick post rwd delivery)
                ae_time = min(medass_cell{15}(medass_cell{15}>sum(trl_mtx(current_trial, [1 11]))));
            case 7 % Trial end
                ae_time = sum(trl_mtx(current_trial, [1 12]));
            otherwise
                error('incorrect alignment input')
        end

        if ~isempty(ae_time) && ~isnan(ae_time) % not all trials have rewards

            % time window bounds
            tw_bounds_trl = [ae_time-tw_bounds(1) ae_time+tw_bounds(2)];

            % neuron activity observed
            frametime_idx = all_pkl_frames>=tw_bounds_trl(1) & all_pkl_frames<=tw_bounds_trl(2) & trl_idx' == itrl;
            frametimes_obs = all_pkl_frames(frametime_idx);
            activity_obs = traces(current_cell, frametime_idx);
            
            % neuron activity interpolated
            frametimes_int = linspace(tw_bounds_trl(1), tw_bounds_trl(2), num_samps); 

            % frametimes_int
            
            if size(frametimes_obs,2) == 0
                continue
            end

            
            activity_int = interp1(frametimes_obs, activity_obs, frametimes_int, 'linear');
            
            % load output
            act_mtx(ic, :, itrl) = activity_int;

        end
        
    end
end
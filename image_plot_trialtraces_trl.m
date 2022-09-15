function image_plot_trialtraces_trl(trl_mtx, all_pkl_frames, trl_idx, traces, neurons, trials, varargin)
% plot each trial as a serperate figure; include all neurons sorted
% according to srt_idx (varargin)

% remove abberations
for i = 1:size(traces)
    traces(i,:) = remove_abberations(traces(i,:), 2.5);
end

% sort neurons
if ~isempty(varargin)
   srt_idx = varargin{1};
else
    srt_idx = 1:length(neurons);
end
neurons = neurons(srt_idx);

% normalize all traces to [0 1]
traces = traces(neurons,:);
traces = traces - min(traces, [], 2);
traces = traces ./ max(traces, [], 2);

% figure
figure; hold on; 

% colors
colors = distinguishable_colors(length(trials));

% iterate through all neurons
for ic = 1:length(neurons)
    
    % iterate through all trials
    for itrl = 1:length(trials)
    
        % current items
        current_trial = trials(itrl);
        
        % 1 cell per row
        %yaxis_row = length(neurons) - ic;
        yaxis_row = length(neurons) - ic + 0.5;
        
        % time
        x_time = all_pkl_frames(trl_idx==current_trial);
        min_time = min(x_time);
        x_time = x_time - min_time; % set start to 0s
        %max_time = max(x_time);
        %x_time = x_time./max_time; % NO
        
        
        % NEURON ACTIVITY
        y_trace = traces(ic, trl_idx==current_trial);
        
        % plot traces
        h = plot(x_time, y_trace.*1.3 + yaxis_row, 'color', .3.*[1 1 1]);
        
        
        
        % plot time in NP
        np_onset = sum(trl_mtx(current_trial, [1 6])) - min_time;
        np_offset = sum(trl_mtx(current_trial, [1 9])) - min_time;
        plot(np_offset.*[1 1], [0 1] + yaxis_row, 'r-', 'linewidth', 2)
        %patch([np_onset np_offset np_offset np_onset], [0 0 1 1] + yaxis_row, 0.7.*[1 1 1], 'FaceAlpha', .3)
        
        % plot head entry
        fd_onset = sum(trl_mtx(current_trial, [1 10])) - min_time;
        fd_offset = sum(trl_mtx(current_trial, [1 10])) - min_time + 2;
        %patch([fd_onset fd_offset fd_offset fd_onset], [0 0 1 1] + yaxis_row, 0.7.*[1 1 1], 'FaceAlpha', .3)
        plot(fd_onset.*[1 1], [0 1] + yaxis_row, 'r-', 'linewidth', 2)
        
        % plot tone on
        tone_on_time = sum(trl_mtx(current_trial, [1 7])) - min_time;
        plot(tone_on_time.*[1 1], [0 1] + yaxis_row, 'r-', 'linewidth', 2)
        
        
        % plot start of random delay
        fd_offset = sum(trl_mtx(current_trial, [1 10])) - min_time + 2;
        plot(fd_offset.*[1 1], [0 1] + yaxis_row, 'r-', 'linewidth', 2)
     
        
        % plot reward time
        rwd_delivery_time = sum(trl_mtx(current_trial, [1 11])) - min_time;
        plot(rwd_delivery_time.*[1 1], [0 1] + yaxis_row, 'r-', 'linewidth', 2)
        
        % plot reward receipt time
        %{
        if ~isnan(sum(trl_mtx(current_trial, [1 11])))
            first_lick_time = min(medass_cell{15}(medass_cell{15}>sum(trl_mtx(current_trial, [1 11])))) - min_time;
            plot(first_lick_time.*[1 1], [0 1] + yaxis_row, '-', 'color', colors(itrl,:), 'linewidth', 2)
        end
        %}
        
        
        
        
        
    end
end
  
% aesthetics
%set(gcf,'Position', [550 55 1010 1283])
set(gca,'TickLength',[0, 0]); box off;
ylim([0.5 length(neurons)+0.5])
ytl_hold = length(neurons):-1:1;
yticklabels(ytl_hold(yticks))
ylabel('Neuron')
xlabel('Time (s)')
xlim([0 max(x_time)])




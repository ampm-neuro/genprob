function image_plot_trialtraces_trl_hm(trl_mtx, all_pkl_frames, trl_idx, traces, neurons, trials, varargin)
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

% iterate through all trials
for itrl = 1:length(trials)
    
        % figure
        figure; hold on; 
    
        % current items
        current_trial = trials(itrl);

        % time
        x_time = all_pkl_frames(trl_idx==current_trial);
        min_time = min(x_time);
        x_time = x_time - min_time; % set start to 0s
        
        % plot activity
        trial_activity = traces(:, trl_idx==current_trial);
        imagesc(trial_activity)
        
        % aesthetics
        set(gca,'TickLength',[0, 0]); box off;
        set(gca, 'YDir','reverse')
        ylim([0.5 length(neurons)+0.5]);
        ytl_hold = length(neurons):-1:1;
        
        %yticklabels(ytl_hold(yticks))
        ylabel('Neuron')
        xlabel('Time (s)')
        xlim([0.5 size(trial_activity,2)])
        title(['Trial ' num2str(current_trial)])
        xticklabels(xticks./100)
        
        
        
        % plot time in NP
        np_offset = sum(trl_mtx(current_trial, [1 9])) - min_time;
        np_offset = np_offset*100 + 0.5;
        plot(np_offset.*[1 1], ylim, 'r-', 'linewidth', 2)
        
        % plot head entry
        fd_onset = sum(trl_mtx(current_trial, [1 10])) - min_time;
        fd_onset = fd_onset*100 + 0.5;
        plot(fd_onset.*[1 1], ylim, 'r-', 'linewidth', 2)
        
        % plot tone on
        tone_on_time = sum(trl_mtx(current_trial, [1 7])) - min_time;
        tone_on_time = tone_on_time*100 + 0.5;
        plot(tone_on_time.*[1 1], ylim, 'r-', 'linewidth', 2)
        
        % plot start of random delay
        fd_offset = sum(trl_mtx(current_trial, [1 10])) - min_time + 2;
        fd_offset = fd_offset*100 + 0.5;
        plot(fd_offset.*[1 1], ylim, 'r-', 'linewidth', 2)
        
        % plot reward time
        rwd_delivery_time = sum(trl_mtx(current_trial, [1 11])) - min_time;
        rwd_delivery_time = rwd_delivery_time*100 + 0.5;
        plot(rwd_delivery_time.*[1 1], ylim, 'r-', 'linewidth', 2)
        
        
        
end

  



%patch([np_onset np_offset np_offset np_onset], [0 0 1 1] + yaxis_row, 0.7.*[1 1 1], 'FaceAlpha', .3)




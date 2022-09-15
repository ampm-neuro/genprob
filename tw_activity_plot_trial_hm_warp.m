function hm_cell = tw_activity_plot_trial_hm_warp(trl_mtx, trl_idx, all_pkl_frames, traces, neurons, trials, tw_bounds, tses)
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

%timewarp traces
[traces, all_pkl_frames, trl_idx, trl_mtx] = ...
    timewarp_traces(traces, all_pkl_frames, trl_idx, trl_mtx, tses);


%activity matrix
act_mtx = tw_activity_II(trl_mtx, trl_idx, all_pkl_frames, traces, neurons, trials, 3, tw_bounds);


% iterate through each neuron
hm_cell = cell(1, length(neurons));
for ic = 1:length(neurons)
    
    % iterate through tones
    for itrl = length(trials):-1:1
        % compute trace
        trace_local = act_mtx(ic, :, itrl);
        trace_local = norm_mtx(trace_local);
        hm_cell{ic}(itrl,:) = trace_local;
    end
    
    % plot
    %figure; 
    hold on
    imagesc(hm_cell{ic})
    
    % aesthetics
    set(gca,'TickLength',[0, 0]); box off;
    set(gca, 'YDir','reverse')
    ylim([0.5 length(trials)+0.5])
    xlim([0 inf])
    title(['Neuron number ' num2str(neurons(ic))])
    
    if rem(length(trials),10) > 1
        yticks([1 rem(length(trials),10):10:length(trials)])
    else
        yticks([1:10:length(trials)])
    end
    yticklabels(fliplr([1 10:10:length(trials)-1 length(trials)]))
    xticks([1 size(act_mtx,2)])
    xticklabels([-tw_bounds(1) tw_bounds(2)])
    xlabel('Time (s)')
    ylabel('Trial')
    colorbar
    
    % red lines
    %
    if ic==1
        align_xval = size(act_mtx,2).*(tw_bounds(1)/sum(tw_bounds)) - 0.5;
        tses = [0 tses];
        tses = tses.*100; % 100 frames per second
        ctses = cumsum(tses);
        ctses = ctses + (align_xval - ctses(4));
    end
    for rl = ctses([1 2 3 4 6]) %EDITED
        plot(rl.*[1 1], ylim, 'r-', 'linewidth', 2)
    end
    plot(ctses(5).*[1 1], ylim, 'r--', 'linewidth', 2)
    hold off
    %}
    
    %figure; 
    %plot(nanmean(flipud(hm_cell)))
    
end

    
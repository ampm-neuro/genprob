function tw_activity_plot_trial_hm(trl_mtx, trl_idx, medass_cell, frame_times, traces, neurons, trials, align_event, tw_bounds)
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

%activity matrix
act_mtx = tw_activity(trl_mtx, trl_idx, medass_cell, frame_times, traces, neurons, trials, align_event, tw_bounds);

% constrain based on input
trl_mtx_local = trl_mtx(trials, :);

% iterate through each neuron
for ic = 1:length(neurons)
    
    % iterate through tones
    trl_mtx = nan(length(trials),size(act_mtx,2));
    for itrl = length(trials):-1:1
    
        % compute trace
        trace_local = act_mtx(ic, :, itrl);
         trace_local = norm_mtx(trace_local);
        trl_mtx(itrl,:) = trace_local;

    end
    
    % plot
    %figure; 
    hold on
    imagesc(flipud(trl_mtx))
    
    % aesthetics
    set(gca,'TickLength',[0, 0]); box off;
    ylim([0.5 length(trials)+0.5])
    xlim([0 inf])
    title(['Neuron number ' num2str(neurons(ic))])
    yticks([1 rem(length(trials),10):10:length(trials)])
    yticklabels(fliplr([1 10:10:length(trials)-1 length(trials)]))
    xticks([1 size(act_mtx,2)])
    xticklabels([-tw_bounds(1) tw_bounds(2)])
    xlabel('Time (s)')
    ylabel('Trial')
    colorbar
    
    % red line
    xval = size(act_mtx,2).*(tw_bounds(1)/sum(tw_bounds)) - 0.5;
    plot(xval.*[1 1], ylim, 'r-')
    
    
    %figure; 
    %plot(nanmean(flipud(trl_mtx)))
end

    
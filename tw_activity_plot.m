function tw_activity_plot(trl_mtx, medass_cell, all_pkl_frames, traces, neurons, trials, align_event, tw_bounds)
% runs tw_activity and produces one plot for each neuron
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
act_mtx = tw_activity(trl_mtx, medass_cell, all_pkl_frames, traces, neurons, trials, align_event, tw_bounds);

% iterate through each neuron
for ic = 1:length(neurons)
    
    % figure
    figure; hold on
    
    % compute mean and se
    activity_mean = nanmean(act_mtx(ic, :, :), 3);
    activity_se = nanstd(act_mtx(ic, :, :), [], 3) ./ sqrt(sum(~isnan(act_mtx),3));
    
    % plot
    xvals = linspace(-tw_bounds(1), tw_bounds(2), length(activity_mean));
    
    plot(xvals, activity_mean, 'linewidth', 2, 'color', 0.0.*[1 1 1])
    plot(xvals, activity_mean+activity_se, 'linewidth', 1, 'color', 0.6.*[1 1 1])
    plot(xvals, activity_mean-activity_se, 'linewidth', 1, 'color', 0.6.*[1 1 1])

    % aesthetics
    set(gca,'TickLength',[0, 0]); box off;
    title(['Neuron number ' num2str(neurons(ic))])
    plot([0 0], ylim, 'r-')
    
end
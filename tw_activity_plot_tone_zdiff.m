function tw_activity_plot_tone_zdiff(trl_mtx, medass_cell, all_pkl_frames, traces, neurons, trials, align_event, tw_bounds)
% runs tw_activity and produces one plot for each neuron comparing activity
% for each of two unique tones
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
% the tw starts 2 seconds before and ends 3s after the align event.
% compute all activity

%activity matrix
act_mtx = tw_activity(trl_mtx, medass_cell, all_pkl_frames, traces, neurons, trials, align_event, tw_bounds);

% constrain based on input
trl_mtx_local = trl_mtx(trials, :);

% unique tones
utones = unique(floor(trl_mtx(:,2)));

% colors
%colors = winter(length(utones));

% legend
legend_input = cell(1,length(utones));
for i = 1:length(utones)
    legend_input{i} = num2str(utones(i));
end

% iterate through each neuron
for ic = 1:length(neurons)
    
    % figure
    if length(neurons)>1
        figure; 
    end
    hold on
    
    % legend
    %legend_trick(colors, '-')
    
    % iterate through both tones
    activity_means = nan(2, size(act_mtx,2));
    activity_stds = nan(2, size(act_mtx,2));
    
    for itone = 1:2
    
        % compute mean and se
        activity_means(itone,:) = nanmean(act_mtx(ic, :, floor(trl_mtx_local(:,2))==utones(itone)), 3);
        activity_stds(itone,:) = nanstd(act_mtx(ic, :, floor(trl_mtx_local(:,2))==utones(itone)), [], 3);

    end

    % compute zdiff
    mean_diffs = abs(activity_means(1,:)-activity_means(2,:));
    mean_stds = mean(activity_stds);
    zdiffs = mean_diffs./mean_stds;
    
    % plot
    xvals = linspace(-tw_bounds(1), tw_bounds(2), size(activity_means,2));
    plot(xvals, zdiffs, 'linewidth', 2)

    % aesthetics
    ylim([0 1])
    set(gca,'TickLength',[0, 0]); box off;
    title(['Neuron number ' num2str(neurons(ic))])
    plot([0 0], ylim, 'r-')
    legend(legend_input);
    
end

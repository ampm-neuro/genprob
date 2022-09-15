function tw_activity_plot_tone(trl_mtx, trl_idx, medass_cell, frame_times, traces, neurons, trials, align_event, tw_bounds)
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
% the tw starts 2 seconds before and ends 3s after the align event.
% compute all activity

%activity matrix
act_mtx = tw_activity_II(trl_mtx, trl_idx, frame_times, traces, neurons, trials, align_event, tw_bounds);

% constrain based on input
trl_mtx_local = trl_mtx(trials, :);

% unique tones
utones = unique(floor(trl_mtx(:,2)));

% colors
colors = winter(length(utones));

% legend
legend_input = cell(1,length(utones));
for i = 1:length(utones)
    legend_input{i} = num2str(utones(i));
end

% iterate through each neuron
for ic = 1:length(neurons)
    
    % figure
    figure; hold on
    legend_trick(colors, '-')
    
    % iterate through tones
    for itone = 1:length(utones)

        %
        % compute mean and se
        activity_mean = nanmean(act_mtx(ic, :, floor(trl_mtx_local(:,2))==utones(itone)), 3);
        activity_se = nanstd(act_mtx(ic, :, floor(trl_mtx_local(:,2))==utones(itone)), [], 3) ./ sqrt(sum(floor(trl_mtx_local(:,2))==utones(itone)));
%
        % plot
        xvals = linspace(-tw_bounds(1), tw_bounds(2), length(activity_mean));

        
        
        %}
        
        %{
        xvals = linspace(-tw_bounds(1), tw_bounds(2), length(activity_mean));
        for itrl = 1:size(act_mtx,3)
            if floor(trl_mtx_local(itrl,2))==utones(1)
                plot(xvals, act_mtx(ic, :, itrl), 'k')
            elseif floor(trl_mtx_local(itrl,2))==utones(2)
                plot(xvals, act_mtx(ic, :, itrl), 'b')
            else
                %error
            end
        end
        %}
        
        
        plot(xvals, activity_mean, 'linewidth', 2, 'color', colors(itone,:))
        plot(xvals, activity_mean+activity_se, 'linewidth', 1, 'color', 0.6.*colors(itone,:))
        plot(xvals, activity_mean-activity_se, 'linewidth', 1, 'color', 0.6.*colors(itone,:))
        
    end
    
    % aesthetics
    ylim([0 1])
    set(gca,'TickLength',[0, 0]); box off;
    title(['Neuron number ' num2str(neurons(ic))])
    plot([0 0], ylim, 'r-')
    legend(legend_input);
    
end

    
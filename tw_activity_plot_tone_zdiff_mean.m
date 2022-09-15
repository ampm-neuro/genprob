function tw_activity_plot_tone_zdiff_mean(trl_mtx, trl_idx, all_pkl_frames, traces, neurons, trials, align_event, tw_bounds, event_time_spacing, varargin)
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

if ~isempty(varargin)
    srt_idx = varargin{1};
else
    srt_idx = 1:length(neurons);
end


%activity matrix
act_mtx = tw_activity_II(trl_mtx, trl_idx, all_pkl_frames, traces, neurons, trials, align_event, tw_bounds);

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

% preallocate
zdiffs = nan(length(neurons), size(act_mtx,2));

% iterate through each neuron
for ic = 1:length(neurons)

    % iterate through both tones
    activity_means = nan(2, size(act_mtx,2));
    activity_stds = nan(2, size(act_mtx,2));
    
    %{
    for itone = 1:2
    
        % compute mean and se
        activity_means(itone,:) = nanmean(act_mtx(ic, :, floor(trl_mtx_local(:,2))==utones(itone)), 3);
        activity_stds(itone,:) = nanstd(act_mtx(ic, :, floor(trl_mtx_local(:,2))==utones(itone)), [], 3);

    end

    % compute zdiff
    mean_diffs = abs(activity_means(1,:)-activity_means(2,:));
    mean_stds = mean(activity_stds);
    %zdiffs(ic,:) = mean_diffs./mean_stds;
    zdiffs(ic,:) = computeCohen_d(;
    
    %}
    
    for itw = 1:size(zdiffs,2)
        
        activity_t1 = act_mtx(ic, itw, floor(trl_mtx_local(:,2))==utones(1));
        activity_t2 = act_mtx(ic, itw, floor(trl_mtx_local(:,2))==utones(2));
        
        activity_t1 = squeeze(activity_t1);
        activity_t2 = squeeze(activity_t2);
        
        if all([sum(~isnan(activity_t1))>=(0.8*length(activity_t1)) sum(~isnan(activity_t2))>=(0.8*length(activity_t2))])
            zdiffs(ic,itw) = abs(computeCohen_d(activity_t1, activity_t2));
        end
    end
    
end

hold on

% mean and se
zdiffs_mean = nanmean(zdiffs);
zdiffs_se = nanstd(zdiffs, [], 1)./sqrt(nansum(zdiffs,1));


% plot
xvals = linspace(-tw_bounds(1), tw_bounds(2), size(activity_means,2));
h1 = plot(xvals, zdiffs_mean, 'linewidth', 2);
plot(xvals, zdiffs_mean+zdiffs_se, 'linewidth', 1, 'color', get(h1,'Color'))
plot(xvals, zdiffs_mean-zdiffs_se, 'linewidth', 1, 'color', get(h1,'Color'))

% aesthetics
ylim([0 1])
set(gca,'TickLength',[0, 0]); box off;
plot([0 0], ylim, 'r-')
legend(legend_input);


%% plot zdiff hm

%rvals_all = abs(rvals_all);
zdiffs = zdiffs(srt_idx,:);

% remove empty columns
zdiffs = zdiffs(:,~isnan(sum(zdiffs)));

figure; imagesc(zdiffs)


% red lines
event_frame = cumsum([0 event_time_spacing]).*100;

xticks_hold = [];
for i = event_frame 
   plot(i.*[1 1], ylim, 'r-', 'linewidth', 1.5) 
   xticks_hold = [xticks_hold i];
end

% aesthetics
set(gca, 'YDir','reverse')
ylim([0+0.5 size(rvals_hm,1)+0.5])
xlim([0.5 size(rvals_hm,2)+0.5])
xticks(xticks_hold)
xticklabels_universe = {'NP1', 'NP2', 'HE1', 'TO', 'RD', 'RWD', 'TE'};
xticklabels(xticklabels_universe)

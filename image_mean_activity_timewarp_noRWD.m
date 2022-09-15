function [all_mtx, sort_idx, sort_max_pix, max_pix, all_mtx_unsort, hm_cells] = image_mean_activity_timewarp_noRWD(trl_mtx, trl_idx, frame_times, traces, neurons, trials, time_series_event_spacing, varargin)
% warps the activity timeline of every trial and then plots the mean trial 
% activity of each cell
%
% rewarded and unrewarded trials plotted differently?
% different for each tone?

% sort index input
if ~isempty(varargin)
   sort_idx_input = varargin{1}; 
end

% fixed trial timeline: nose_poke onset (0, implied), nose_poke offset, head_entry, tone
% on, trial start, reward, trial end
%time_series_event_spacing = [0.2 1.1 1.5 2.0 3.5 2.0];

% warp time
[warp_traces, warp_frame_times, warp_trl_idx, warp_trl_mtx] = ...
    timewarp_traces(traces, frame_times, trl_idx, trl_mtx, time_series_event_spacing);

% compute cell activity on each trial
hm_cells = tw_activity_trial_hm_full_warp(warp_trl_mtx, warp_trl_idx,...
    warp_frame_times, warp_traces, neurons, trials,...
    [sum(time_series_event_spacing(1:2))+5 sum(time_series_event_spacing(3:end))+5], time_series_event_spacing);

% display number of trials
if size(hm_cells{1},1)<2
    disp(['TOO FEW TRIALS: ' num2str(size(hm_cells{1},1)) ' trials'])
    all_mtx = [];
    sort_idx = [];
    return
else
    disp([num2str(size(hm_cells{1},1)) ' trials'])
end

% average activity
all_mtx = []; 
for ihm = 1:length(hm_cells) 
    hm_cells{ihm} = hm_cells{ihm}(:,sum(isnan(hm_cells{ihm}))==0);
    all_mtx = [all_mtx; nanmean(hm_cells{ihm})]; 
end

% remove empty columns
all_mtx = all_mtx(:,nansum(all_mtx)>0);

% chop off edges?
edge_chop = 0;
%all_mtx = all_mtx(:, edge_chop+1:end-edge_chop);
all_mtx(:, [1:edge_chop+1 end-edge_chop:end]) = nan;

% smooth rows
for i = 1:2
    for irow = 1:size(all_mtx,1)
        all_mtx(irow,:) = smooth(all_mtx(irow,:),5);
    end
end

% normalize rows
all_mtx = (all_mtx-nanmin(all_mtx,[],2))./nanmax((all_mtx-min(all_mtx,[],2)),[],2);

% find max pixel
%
max_pix = nan(size(all_mtx,1),1); 
for irow = 1:size(all_mtx,1)
    max_pix(irow) = find(all_mtx(irow,:)==max(all_mtx(irow,:))); 
end
[sort_max_pix, sort_idx] = sort(max_pix);

% sort rows
all_mtx_unsort = all_mtx;
if exist('sort_idx_input', 'var')
    all_mtx = all_mtx(sort_idx_input,:);
    sort_max_pix = max_pix(sort_idx_input);
else
    all_mtx = all_mtx(sort_idx,:);
end
%}

% plot
figure; hold on
imagesc(all_mtx)

% red line
event_frame = cumsum([0 time_series_event_spacing]).*100;

xticks_hold = [];
for i = event_frame 
   plot(i.*[1 1], ylim, 'r-', 'linewidth', 1.5) 
   xticks_hold = [xticks_hold i];
end

% aesthetics
set(gca, 'YDir','reverse')
ylim([0+0.5 size(all_mtx,1)+0.5])
xlim([0.5 size(all_mtx,2)+0.5])
xticks(xticks_hold)
xticklabels_universe = {'NP1', 'NP2', 'HE1', 'TO', 'RD', 'RWD', 'TE'};
xticklabels(xticklabels_universe)

% peaks
for i = 1:size(all_mtx,1)
    plot(sort_max_pix(i), i, 'ko')
end

% lineplot
figure; hold on 
plot(smooth(mean(all_mtx), 10), 'k-', 'linewidth', 2)
ylim([0 1])
set(gca,'TickLength',[0, 0]); box off;

% red line
event_frame = cumsum([0 time_series_event_spacing]).*100;
for i = event_frame 
   plot(i.*[1 1], ylim, 'r-', 'linewidth', 1.5) 
end
xlim([0.5 size(all_mtx,2)+0.5])





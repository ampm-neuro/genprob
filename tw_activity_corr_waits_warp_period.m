function [rvals_hm] = tw_activity_corr_waits_warp_period(trl_mtx, trl_idx, frame_times, traces, neurons, align_event, tw_bounds, event_time_spacing, varargin)
% computes correlation between activity in each task period and wait times
% during probe trials

if ~isempty(varargin)
    srt_idx = varargin{1};
else
    srt_idx = 1:length(neurons);
end

% warp
[warp_traces, warp_frame_times, warp_trl_idx, warp_trl_mtx] = timewarp_traces(traces, frame_times, trl_idx, trl_mtx, event_time_spacing);

%activity matrix
trials = unique(trl_idx);
act_mtx = tw_activity_II(warp_trl_mtx, warp_trl_idx, warp_frame_times, warp_traces, neurons, trials, align_event, tw_bounds);


% REMOVE TRIALS WITHOUT VIDEO
trl_mtx = trl_mtx(unique(trl_idx),:);

% smooth activity mtx
smwdw = 15;
for icell = 1:size(act_mtx,1)
    for itrl = 1:size(act_mtx,3)
        nnan_idx = ~isnan(act_mtx(icell,:,itrl));
        act_mtx(icell, nnan_idx, itrl) = smooth(act_mtx(icell, nnan_idx, itrl), smwdw);
    end
end

% average across window
act_mtx_mean = nanmean(act_mtx,2);
act_mtx_mean = reshape(act_mtx_mean, size(act_mtx,1), size(act_mtx,3));

% probe trials only
probe_trial_nums = find(trl_mtx(:,3)==0);
trl_mtx_p = trl_mtx(probe_trial_nums,:);
act_mtx_mean = act_mtx_mean(:,probe_trial_nums);
act_mtx = act_mtx(:,:,probe_trial_nums);
probe_waits = trl_mtx_p(:,12);

% compute r, p , and plot fit line for each cell
%{
rpvals = nan(size(act_mtx,1), 2); 
for ic = 1:size(act_mtx_mean,1)
    
    % plot fit line
    figure; 
    [r, p] = fit_line(act_mtx_mean(ic,:)', probe_waits); 
    title(['neuron number ' num2str(ic) '; r=' num2str(r) ', p=' num2str(p)]);
    
    % aesthetics
    xlabel('mean activity over window')
    ylabel('wait time')
    
    % load values
    rpvals(ic,:) = [r p]; 
     
end
%}


% compute average rvalue at each time point over window
%

% preallocate rval matrix
rvals_all = nan(size(act_mtx,1), size(act_mtx,2));

% over all cells
for ic = 1:size(act_mtx,1)
    
    % over all time windows
    for iw = 1:size(act_mtx,2)
        % load values
        rvals_all(ic,iw) = corr(squeeze(act_mtx(ic,iw,:)), probe_waits); 
    
    end
    
end

%rvals_all = abs(rvals_all);
rvals_hm = rvals_all(srt_idx,:);

% remove empty columns
rvals_hm = rvals_hm(:,~isnan(sum(rvals_hm)));

% smooth rows
for i = 1:2
    for irow = 1:size(rvals_hm,1)
        rvals_hm(irow,:) = smooth(rvals_hm(irow,:),5);
    end
end


% heatmap
%
figure; hold on;
imagesc(rvals_hm)
%}

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



% average
figure; 
hold on;
h = plot(nanfastsmooth(nanmean(rvals_all),11), 'linewidth', 2);
plot(nanfastsmooth(nanmean(rvals_all)+nanstd(rvals_all)./sqrt(size(rvals_all,1)),11), '-', 'color', get(h,'Color') , 'linewidth', 1)
plot(nanfastsmooth(nanmean(rvals_all)-nanstd(rvals_all)./sqrt(size(rvals_all,1)),11), '-', 'color', get(h,'Color') , 'linewidth', 1)
plot(tw_bounds(1).*100.*[1 1], ylim, 'r-', 'linewidth', 2)




    
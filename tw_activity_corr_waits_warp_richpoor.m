function [rvals_hm] = tw_activity_corr_waits_warp_richpoor(trl_mtx, trl_idx, frame_times, traces, neurons, align_event, tw_bounds, event_time_spacing, varargin)
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

if length(varargin)==1
    srt_idx = varargin{1};
elseif length(varargin)==2
    srt_idx = varargin{1};
    peaks = varargin{2};
else
    srt_idx = 1:length(neurons);
    peaks = nan(size(srt_idx));
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

% probe trials only
probe_trial_nums = find(trl_mtx(:,3)==0);
trl_mtx_p = trl_mtx(probe_trial_nums,:);
act_mtx = act_mtx(:,:,probe_trial_nums);
probe_waits = trl_mtx_p(:,12);


% compute average rvalue at each time point over window for rich and then
% poor probes
%
itrl_set{1} = find(rich_trl_idx(trl_mtx_p));
itrl_set{2} = find(~rich_trl_idx(trl_mtx_p));
for iis = 1:2

    % preallocate rval matrix
    rvals_all = nan(size(act_mtx,1), size(act_mtx,2));

    % over all cells
    for ic = 1:size(act_mtx,1)

        % over all time windows
        for iw = 1:size(act_mtx,2)
            % load values
            rvals_all(ic,iw) = corr(squeeze(act_mtx(ic,iw,itrl_set{iis})), probe_waits(itrl_set{iis})); 

        end

    end

    %rvals_all = abs(rvals_all);
    rvals_hm = rvals_all(srt_idx,:);

    % remove empty columns
    rvals_hm = rvals_hm(:,~isnan(sum(rvals_hm)));

    % smooth rows
    for i = 1:2
        for irow = 1:size(rvals_hm,1)
            rvals_hm(irow,:) = smooth(rvals_hm(irow,:),10);
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

    % peaks plot
    for i = 1:size(rvals_hm,1)
        plot(peaks(i), i, 'ko')
    end
    
    % aesthetics
    set(gca, 'YDir','reverse')
    ylim([0+0.5 size(rvals_hm,1)+0.5])
    xlim([0.5 size(rvals_hm,2)+0.5])
    xticks(xticks_hold)
    xticklabels_universe = {'NP1', 'NP2', 'HE1', 'TO', 'RD', 'RWD', 'TE'};
    xticklabels(xticklabels_universe)
    if iis==1; title('rich probes'); else; title('poor probes'); end

    
    
    rvals_all = rvals_hm;
    
    % average
    figure; 
    hold on;
    h = plot(nanfastsmooth(nanmean(rvals_all),11), 'linewidth', 2);
    plot(nanfastsmooth(nanmean(rvals_all)+nanstd(rvals_all)./sqrt(size(rvals_all,1)),11), '-', 'color', get(h,'Color') , 'linewidth', 1)
    plot(nanfastsmooth(nanmean(rvals_all)-nanstd(rvals_all)./sqrt(size(rvals_all,1)),11), '-', 'color', get(h,'Color') , 'linewidth', 1)
    plot(tw_bounds(1).*100.*[1 1], ylim, 'r-', 'linewidth', 2)
    xticks(xticks_hold)
    xticklabels_universe = {'NP1', 'NP2', 'HE1', 'TO', 'RD', 'RWD', 'TE'};
    xticklabels(xticklabels_universe)
    if iis==1; title('rich probes'); else; title('poor probes'); end
    
    
    % corr around firing rate peak
    %
    tw = 2;
    peak_FR_corr = nan(size(rvals_all,1),tw*100+1);
    bins_out_meta = 100*(tw/2);
    for ir = 1:size(rvals_all,1)
        bins_out_lo = bins_out_meta;
        bins_out_hi = bins_out_meta;

        % avoid edges
        if peaks(ir)==0
            continue
        elseif peaks(ir) - bins_out_lo < 0
            bins_out_lo = peaks(ir)-1;
        end
        if peaks(ir) >= size(rvals_all,2)
            continue
        elseif peaks(ir)+bins_out_hi > size(rvals_all,2)
            bins_out_hi = (size(rvals_all,2)-peaks(ir))-1;
        end

        % load
        peak_FR_corr(ir,((bins_out_meta + 1)-bins_out_lo):((bins_out_meta + 1)+bins_out_hi)) = rvals_all(ir, peaks(ir)-bins_out_lo : peaks(ir)+bins_out_hi);
    end
    peak_r_mean = nanmean(peak_FR_corr);
    peak_r_ste = nanstd(peak_FR_corr)./sqrt(sum(~isnan(peak_FR_corr)));
    
    
    figure; imagesc(peak_FR_corr)
    
    figure; hold on;
    plot(1:size(peak_FR_corr,2), peak_r_mean, 'k-', 'linewidth', 2)
    plot(1:size(peak_FR_corr,2), peak_r_mean-peak_r_ste, 'k-', 'linewidth', 1)
    plot(1:size(peak_FR_corr,2), peak_r_mean+peak_r_ste, 'k-', 'linewidth', 1)
    
    xlim([0 bins_out_meta*2+2])
    xticks([1 (bins_out_meta*2+2)/2 bins_out_meta*2+1])
    xticklabels({num2str(-tw/2), 'Mean rate peak', num2str(tw/2)})
    ylim([-1 1])
    plot(xlim, [0 0], 'k--')
    plot((bins_out_meta*2+2)/2.*[1 1], ylim, 'r-')
    
   %}

   
   
end




    
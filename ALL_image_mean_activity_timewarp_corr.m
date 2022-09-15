function RD_corr_sum = ALL_image_mean_activity_timewarp_corr(file_keywords, first_last)
% combine cells from a single stage from all animals

%% prep
% fixed trial timeline: nose_poke onset (0, implied), nose_poke offset, head_entry, tone
% on, start random delay, reward or abandon, trial end (if rewarded)
time_series_event_spacing = [0.2 1.1 1.0 2.0 8.5 2.0];

% get all files
fp = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_imaging_hpc\';
session_files = get_file_paths_targeted(fp, file_keywords);

% find unique subjects
subjs = [];
for isf = 1:size(session_files,1)
    subjs = [subjs; {session_files{isf}(strfind(session_files{isf}, '\\')+2 : strfind(session_files{isf}, '\\')+9)}];
end
subjs = unique(subjs);

% find one file per subject
session_files_hold = [];
for isubj = 1:length(subjs)
    subj_sf = session_files(contains(session_files, subjs{isubj}));
    if strcmp(first_last, 'first')
        session_files_hold = [session_files_hold; subj_sf(1)];
    elseif strcmp(first_last, 'last')
        session_files_hold = [session_files_hold; subj_sf(end)];
    elseif strcmp(first_last, 'all')
        session_files_hold = [session_files_hold; subj_sf(:)];
    else
        error('do better with first_last input')
    end
end
session_files_hold
session_files_hold = session_files_hold(contains(session_files_hold, 'LED')); %LED ONLY
%session_files_hold = session_files_hold(contains(session_files_hold, 'preprobe_02'));
%{
session_files_hold = session_files_hold(~contains(session_files_hold, 'preprobe_01'));
%session_files_hold = session_files_hold(~contains(session_files_hold, 'preprobe_02'));
session_files_hold = session_files_hold(~contains(session_files_hold, 'preprobe_03'));
session_files_hold = session_files_hold(~contains(session_files_hold, 'preprobe_04'));
%session_files_hold = session_files_hold(~contains(session_files_hold, 'preprobe_05'));
%session_files_hold = session_files_hold(~contains(session_files_hold, 'preprobe_06'));
%session_files_hold = session_files_hold(~contains(session_files_hold, 'postprobe_04'));
%}

%session_files_hold = session_files_hold(~contains(session_files_hold, '-01.mat'));
session_files = session_files_hold


%% compute
% iterate through session files
all_mtx = [];
all_sort = [];
all_corr_hm = nan(1,sum(time_series_event_spacing(1:end-1).*100));
all_smp = [];
all_mean_act = [];
all_probe_waits = [];
all_rich_idx = [];
for isf = 1:size(session_files,1)

    % load session file
    session_files{isf}
    load(session_files{isf}, 'frame_times', 'traces', 'trl_idx', 'trl_mtx')
    
    %plot_trials_trlmtx_altcolor(trl_mtx)
    
    if ~exist('trl_idx', 'var')
        continue
    end

    % compute matrices and order    
    [sf_mtx,sf_srt_idx,sort_max_pix] = image_mean_activity_timewarp...
        (trl_mtx, trl_idx, frame_times, traces, 1:size(traces,1),...
        intersect(unique(trl_idx), find(trl_mtx(:,3)==0)), time_series_event_spacing);
    
    % close figures
    close; close

    % correlation
    %{
    [corr_hm, mean_activity_trl, probe_waits, rich_trl_idx] = tw_activity_corr_waits_warp(trl_mtx, trl_idx, frame_times, traces,...
        1:size(traces,1), 4, [sum(time_series_event_spacing(1:2)) sum(time_series_event_spacing(3:end))],...
        time_series_event_spacing, sf_srt_idx);
    %}
    [corr_hm, mean_activity_trl, probe_waits, rich_trl_idx] = tw_activity_corr_waits_warp(trl_mtx, trl_idx, frame_times, traces,...
        1:size(traces,1), 4, [sum(time_series_event_spacing(1:2)) sum(time_series_event_spacing(3:end))],...
        time_series_event_spacing);

    % close figures
    close; close

    % load
    all_mtx = [all_mtx; sf_mtx];
    all_sort = [all_sort; sf_srt_idx];
    
    % adjust for small timing errors
    length_dif = size(all_corr_hm,2)-size(corr_hm,2);
    if length_dif > 0 && length_dif < 10
        corr_hm = [corr_hm repmat(corr_hm(:,end), 1, length_dif)];
    elseif length_dif > 5
        error('matrix dimensions wrong ampm')
    end
    
    all_corr_hm = [all_corr_hm; corr_hm(sf_srt_idx,:)]; %sort by activity peak timing
    all_smp = [all_smp; sort_max_pix];
    all_mean_act = [all_mean_act; zscore_mtx(mean_activity_trl)];
    all_probe_waits = [all_probe_waits; zscore_mtx(probe_waits)];
    all_rich_idx = [all_rich_idx; rich_trl_idx];
    
    % clear lingering variables
    clearvars('frame_times', 'traces', 'trl_idx', 'trl_mtx')
end
all_rich_idx = logical(all_rich_idx);


%% plot trial mean activity (during random delay)
figure; hold on
legend_trick([0 0 0], 'k-')
plot(all_mean_act(~all_rich_idx), all_probe_waits(~all_rich_idx), 'ko')
plot(mean(all_mean_act(~all_rich_idx)), mean(all_probe_waits(~all_rich_idx)), 'k.', 'markersize', 50)
plot(all_mean_act(all_rich_idx), all_probe_waits(all_rich_idx), 'ro')
plot(mean(all_mean_act(all_rich_idx)), mean(all_probe_waits(all_rich_idx)), 'r.', 'markersize', 50)
[r,p] = fit_line(all_mean_act, all_probe_waits);
axis square
axis([-5 5 -5 5])
hold on; plot(xlim, [0 0], 'k--'); plot([0 0], ylim, 'k--')
ylabel('probe waits (z)')
xlabel('mean activity (z)')
legend({['r=' num2str(round(r*1000)/1000) ';p=' num2str(round(p*1000)/1000)]}, 'location', 'northwest')
title(file_keywords)

figure; 
errorbar_barplot([{all_mean_act(all_rich_idx)}, {all_mean_act(~all_rich_idx)}])
xticks([1 2])
xticklabels({'rich', 'poor'})
[~, rich_p, ~, rich_stats] = ttest(all_mean_act(all_rich_idx))
title(['activity; p = ' num2str(rich_p)])



%% sort overall
% sort activity peak rows by max pixel
max_pix = nan(size(all_mtx,1),1); 
for irow = 1:size(all_mtx,1)
    max_pix(irow) = find(all_mtx(irow,:)==max(all_mtx(irow,:))); 
end
[~, sort_idx_activity] = sort(max_pix);
all_mtx = all_mtx(sort_idx_activity,:);

% sort correlation rows
all_corr_hm = all_corr_hm(sort_idx_activity,:);
all_smp = all_smp(sort_idx_activity,:);


%% plot activity peak timing
% plot
figure; hold on; title(file_keywords)
imagesc(all_mtx)
ylabel('Neurons')
xlabel('Time')

% red line
event_frame = cumsum([0 time_series_event_spacing]).*100;
xticks_hold = [];
for i = event_frame 
   plot(i.*[1 1], ylim, 'r-', 'linewidth', 1.5) 
   xticks_hold = [xticks_hold i];
end


% black line
plot(all_smp,1:size(all_mtx,1), 'k-', 'linewidth', 2)

% aesthetics
set(gca, 'YDir','reverse')
ylim([0+0.5 size(all_mtx,1)+0.5])
xlim([0.5 size(all_mtx,2)+0.5])
xticks(xticks_hold)
xticklabels_universe = {'NP1', 'NP2', 'HE1', 'TO', 'RD', 'RWD', 'TE'};
xticklabels(xticklabels_universe)

% save
var_name = ['act_peak_timing_hm_' file_keywords{:} '.pdf'] 
print(['C:\Users\ampm1\Desktop\paul_meet\next\' var_name], '-dpdf', '-painters', '-bestfit')

% lineplot mean FR
figure; hold on; title(file_keywords)
plot(smooth(mean(all_mtx), 10), 'k-', 'linewidth', 2)
plot(smooth(mean(all_mtx)+std(all_mtx)./sqrt(size(all_mtx,1)), 10), 'k-', 'linewidth', 1)
plot(smooth(mean(all_mtx)-std(all_mtx)./sqrt(size(all_mtx,1)), 10), 'k-', 'linewidth', 1)
ylim([0 1])
set(gca,'TickLength',[0, 0]); box off;

% red line
event_frame = cumsum([0 time_series_event_spacing]).*100;
for i = event_frame 
   plot(i.*[1 1], ylim, 'r-', 'linewidth', 1.5) 
end
xlim([0.5 size(all_mtx,2)+0.5])
xlabel('Time')

% save
var_name = ['act_peak_timing_line_' file_keywords{:} '.pdf']; 
print(['C:\Users\ampm1\Desktop\paul_meet\next\' var_name], '-dpdf', '-painters', '-bestfit')

%% plot correlation
% plot
figure; hold on; title(file_keywords)
imagesc(all_corr_hm)
caxis([-.8 .8])
colorbar
ylabel('Neurons')
xlabel('Time')

RD_corr_sum = mean(all_corr_hm(:,sum(time_series_event_spacing(1:4))*100:end),2);
RD_corr_se = std(all_corr_hm(:,sum(time_series_event_spacing(1:4))*100:end),[],2);
    RD_corr_se = RD_corr_se./sqrt(size(all_corr_hm,2));
    
% red line
event_frame = cumsum([0 time_series_event_spacing]).*100;
xticks_hold = [];
for i = event_frame 
   plot(i.*[1 1], ylim, 'r-', 'linewidth', 1.5) 
   xticks_hold = [xticks_hold i];
end

% black line
plot(all_smp,1:size(all_corr_hm,1), 'k-', 'linewidth', 2)

% aesthetics
set(gca, 'YDir','reverse')
ylim([0+0.5 size(all_corr_hm,1)+0.5])
xlim([0.5 size(all_corr_hm,2)+0.5])
xticks(xticks_hold)
xticklabels_universe = {'NP1', 'NP2', 'HE1', 'TO', 'RD', 'RWD', 'TE'};
xticklabels(xticklabels_universe)

% save
var_name = ['correlation_hm_' file_keywords{:} '.pdf']; 
print(['C:\Users\ampm1\Desktop\paul_meet\next\' var_name], '-dpdf', '-painters', '-bestfit')

% lineplot
figure; hold on; title(file_keywords)
plot(nanmean(all_corr_hm), 'k-', 'linewidth', 2)
plot(nanmean(all_corr_hm)+nanstd(all_corr_hm)./sqrt(sum(~isnan(all_corr_hm))), 'k-', 'linewidth', 1)
plot(nanmean(all_corr_hm)-nanstd(all_corr_hm)./sqrt(sum(~isnan(all_corr_hm))), 'k-', 'linewidth', 1)
%plot(smooth(nanmean(all_corr_hm)+nanstd(all_corr_hm), 10), 'k-', 'linewidth', 1)
%plot(smooth(nanmean(all_corr_hm)-nanstd(all_corr_hm), 10), 'k-', 'linewidth', 1)
ylim([-1 1])
set(gca,'TickLength',[0, 0]); box off;
xlabel('Time')
hold on; plot(xlim, [0 0], 'k--')
title('mean corr')

% red line
event_frame = cumsum([0 time_series_event_spacing]).*100;
for i = event_frame 
   plot(i.*[1 1], ylim, 'r-', 'linewidth', 1.5) 
end
xlim([0.5 size(all_mtx,2)+0.5])

% save
var_name = ['correlation_line_' file_keywords{:} '.pdf']; 
print(['C:\Users\ampm1\Desktop\paul_meet\next\' var_name], '-dpdf', '-painters', '-bestfit')


% bar plot showing which cells have high r vals DURING random delay
figure; hold on; title(file_keywords)
bar(1:length(RD_corr_sum), RD_corr_sum)
ylim([-0.8 0.8])
xlim([1 length(RD_corr_sum)])
set(gca,'TickLength',[0, 0]); box off;
xlabel('Neurons')
ylabel('rval')

% barplot summary
figure; violin(RD_corr_sum);

% smoothed line version of bar plot above
figure; hold on; title(file_keywords)
[wdw_mean, ~, wdw_se] = nansmooth_ampm(RD_corr_sum,15);
plot(1:length(wdw_mean), wdw_mean, 'k-')
plot(1:length(wdw_mean), wdw_mean+wdw_se, '-', 'color', .7.*[1 1 1])
plot(1:length(wdw_mean), wdw_mean-wdw_se, '-', 'color', .7.*[1 1 1])
ylim([-0.1 0.25])
set(gca,'TickLength',[0, 0]); box off;
xlabel('Neurons')
ylabel('rval')
xlim([1 length(RD_corr_sum)])
hold on; plot(xlim, [0 0], 'k--')

% save
var_name = ['correlation_RD_bar_' file_keywords{:} '.pdf']; 
print(['C:\Users\ampm1\Desktop\paul_meet\next\' var_name], '-dpdf', '-painters', '-bestfit')


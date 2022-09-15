function [common_cell_corrs, common_cell_corrs_means, merge_mtx, common_cell_matrices, common_cell_corrs_cellids] = cell_turnover_timewarp_trials(session_cell, cell_registration, time_series_event_spacing)
% plots average timewarped trials for every registered cell in every
% session, with inactive cells receiving zeros for activity

session_cell

%% compute session matrix for each session
%
session_mtx_cell = cell(size(session_cell));
session_max_cell = cell(size(session_cell));
session_min_cell = cell(size(session_cell));
for isesh = 1:length(session_cell)
   
    % load session
    clearvars traces frame_times 
    session_cell{isesh}
    load(session_cell{isesh})
    if ~exist('traces', 'var')
        continue
    elseif strcmp(session_cell{isesh}, 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_imaging_hpc\690330m2\LED_gen10_mevar04-07.mat')
        continue
    end
    
    % session matrix
    unique_trials = unique(trl_idx);
    if length(unique_trials > 41)
        first_41_trials = unique_trials(1:41);
    end
    [sorted_session_mtx, srt_idx, ~, max_time_bin_unsorted, ~, min_time_bin_unsorted] = image_mean_activity_timewarp(...
    trl_mtx, trl_idx, frame_times, traces, 1:size(traces,1), first_41_trials, time_series_event_spacing);

    % unsort
    unsort_idx = nan(size(srt_idx));
    unsort_idx(srt_idx) = 1:length(srt_idx);
    session_mtx_cell{isesh} = sorted_session_mtx(unsort_idx,:);
    session_max_cell{isesh} = max_time_bin_unsorted;
    session_min_cell{isesh} = min_time_bin_unsorted;
    
end
%save('cttt_683472m3.mat', 'session_mtx_cell', 'session_max_cell', 'session_min_cell')
%}
%load('cttt_683472m3.mat', 'session_mtx_cell', 'session_max_cell', 'session_min_cell')
%disp('... loading cttt_683472m3.mat ...')

%% merge first (or last) appearance of every cell into single matrix
merge_mtx = nan(size(cell_registration,1), size(session_mtx_cell{1},2));
merge_mtx_lo = 1;
merge_mtx_hi = size(cell_registration,1);
max_time_bins = nan(size(cell_registration,1),1);
min_time_bins = nan(size(cell_registration,1),1);
%overall_cell_id = 

% sort by first or last appearance (first = 1, last = 2)
first_last = 2;


% sort by first peak timing during first appearance
if first_last == 1
    disp('merged matrix show first appearance of each cell')
    % for each session
    reliability_corm = [];
    for isesh = 1:size(cell_registration,2)

        if isempty(session_mtx_cell{isesh})
            continue
        end
        
        % unique cells
        if isesh==1
            uc_idx = cell_registration(:,isesh)>0;
            unq_cells = cell_registration(uc_idx, isesh);
        else
            uc_idx = sum(cell_registration(:,1:isesh-1),2)==0 & cell_registration(:,isesh)>0;
            unq_cells = cell_registration(uc_idx, isesh);
        end

        
        % load each unique cell
        merge_mtx(uc_idx, :) = session_mtx_cell{isesh}(unq_cells,:);
        max_time_bins(uc_idx) = session_max_cell{isesh}(unq_cells);
        min_time_bins(uc_idx) = session_min_cell{isesh}(unq_cells);

        % update index
        merge_mtx_lo = merge_mtx_lo+length(unq_cells);
       
        
        %{
        session_mtx_cell_hold{isesh} = session_mtx_cell{isesh}(unq_cells,:);
       
        for icell = 1:size(session_mtx_cell_hold{isesh},1)
            [~, trial_mtx]  = one_cell_reliability(icell, session_mtx_cell_hold);
            reliability_corm = cat(3,reliability_corm, nanmean(trial_mtx,3));
        end
        %}
        
    end
    %reliability_corm = nanmean(reliability_corm,3);
    
    
    
% sort by last peak timing during last appearance
elseif first_last == 2
    disp('merged matrix show last appearance of each cell')
    for isesh = 1:size(cell_registration,2)%size(cell_registration,2):-1:1
        
        if isempty(session_mtx_cell{isesh})
            continue
        end
        
        % unique cells
        if isesh==size(cell_registration,2)
            uc_idx = cell_registration(:,isesh)>0;
            unq_cells = cell_registration(uc_idx, isesh);
        else
            uc_idx = sum(cell_registration(:, isesh+1:1:size(cell_registration,2)),2)==0 & cell_registration(:,isesh)>0;
            unq_cells = cell_registration(uc_idx, isesh);
        end
        
        % load each unique cell
        merge_mtx(uc_idx, :) = session_mtx_cell{isesh}(unq_cells,:);
        max_time_bins(uc_idx) = session_max_cell{isesh}(unq_cells);
        min_time_bins(uc_idx) = session_min_cell{isesh}(unq_cells);

        % update index
        merge_mtx_hi = merge_mtx_hi-length(unq_cells);
    
    end         
end



%% sort according to merged matrix
min_max = 2;
if min_max == 1 % sort by min
    disp('merged matrix sorts by time window with minimum firing rate')
    [max_time_bins_sorted, merge_sort_idx] = sort(min_time_bins);
    cell_registration_sort = cell_registration(merge_sort_idx,:);
    merge_mtx_sorted = merge_mtx(merge_sort_idx,:);
elseif min_max == 2 % sort by max
    disp('merged matrix sorts by time window with maximum firing rate')
    [max_time_bins_sorted, merge_sort_idx] = sort(max_time_bins);
    cell_registration_sort = cell_registration(merge_sort_idx,:);
    merge_mtx_sorted = merge_mtx(merge_sort_idx,:);
end


%% plot merged matrix max
figure; hold on

% plot
imagesc(merge_mtx_sorted); 

% red line
event_frame = cumsum([0 time_series_event_spacing]).*100;
xticks_hold = [];
for i = event_frame 
   plot(i.*[1 1], ylim, 'r-', 'linewidth', 1.0) 
   xticks_hold = [xticks_hold i];
end

% aesthetics
set(gca, 'YDir','reverse')
ylim([0+0.5 size(merge_mtx_sorted,1)+0.5])
xlim([0.5 size(merge_mtx_sorted,2)+0.5])
xticks(xticks_hold)
xticklabels_universe = {'NP1', 'NP2', 'HE1', 'TO', 'RD', 'RWD', 'TE'};
xticklabels(xticklabels_universe)
title('Merged matrix')

% mean activity over time
figure; hold on
plot(nanmean(merge_mtx_sorted))
event_frame = cumsum([0 time_series_event_spacing]).*100;
xticks_hold = [];
for i = event_frame 
   plot(i.*[1 1], ylim, 'r-', 'linewidth', 1.0) 
   xticks_hold = [xticks_hold i];
end

% proportion of peaks
figure; hold on
binsize = 20;
mtbs_hist = max_time_bins_sorted(~isnan(max_time_bins_sorted));
histogram(mtbs_hist, 10:binsize:size(merge_mtx_sorted,2)+realmin, 'normalization', 'probability')
ylim([0 1])
xlim([0 size(merge_mtx_sorted,2)])
event_frame = cumsum([0 time_series_event_spacing]).*100;
xticks_hold = [];
for i = event_frame 
   plot(i.*[1 1], ylim, 'r-', 'linewidth', 1.0) 
   xticks_hold = [xticks_hold i];
end
xticks(xticks_hold)
xticklabels_universe = {'NP1', 'NP2', 'HE1', 'TO', 'RD', 'RWD', 'TE'};
xticklabels(xticklabels_universe)



%% plot average peak and surround
figure; hold on
pas = nan(size(merge_mtx_sorted,1), 251);
for icell = 1:size(merge_mtx_sorted,1)
    if max_time_bins_sorted(icell) >= 20 && max_time_bins_sorted(icell) <= 430
        surround_window = max_time_bins_sorted(icell)-125:max_time_bins_sorted(icell)+125;
        sw_in = 1:length(surround_window);
        sw_in(surround_window<1 | surround_window>size(merge_mtx_sorted,2)) = [];
        surround_window(surround_window<1 | surround_window>size(merge_mtx_sorted,2)) = [];
        pas(icell,sw_in) = merge_mtx_sorted(icell, surround_window);
        %plot(pas(icell,:), '-', 'linewidth', 1, 'color', 0.8.*[1 1 1])
    end
end
plot(nanmean(pas), 'k-', 'linewidth', 3)
plot(nanmean(pas)+nanstd(pas), 'k-', 'linewidth', 1) %./sqrt(sum(~isnan(pas)))
plot(nanmean(pas)-nanstd(pas), 'k-', 'linewidth', 1)
xlim([1 size(pas,2)])
ylim([0 1.01])
hold  on; plot(126.*[1 1], ylim, 'k--')


drawnow
%error('pls stop')

%% plot sorted session matrices

% iterate through sessions
full_session_mtx_3d = nan(size(cell_registration_sort,1), size(session_mtx_cell{isesh},2), length(session_mtx_cell));
figure; 
for isesh = 1:length(session_mtx_cell)
    
    % full session matrix
    full_session_mtx = zeros(size(cell_registration_sort,1), size(session_mtx_cell{isesh},2));
    active_cells_idx = cell_registration_sort(:,isesh)>0;
    if sum(active_cells_idx)==0 || isempty(session_mtx_cell{isesh}) 
        continue; 
    end
    active_cell_session_ids = cell_registration_sort(active_cells_idx, isesh);
    full_session_mtx(active_cells_idx, :) = session_mtx_cell{isesh}(active_cell_session_ids, :);
    full_session_mtx_3d(:,:,isesh) = full_session_mtx;
    
    % plot
    subplot(1,length(session_mtx_cell), isesh); hold on
    imagesc(full_session_mtx); title(num2str(isesh))
    
    % red line
    event_frame = cumsum([0 time_series_event_spacing]).*100;
    xticks_hold = [];
    for i = event_frame 
       plot(i.*[1 1], ylim, 'r-', 'linewidth', 1.0) 
       xticks_hold = [xticks_hold i];
    end

    % aesthetics
    set(gca, 'YDir','reverse')
    ylim([0+0.5 size(full_session_mtx,1)+0.5])
    xlim([0.5 size(full_session_mtx,2)+0.5])
    xticks(xticks_hold)
    xticklabels_universe = {'NP1', 'NP2', 'HE1', 'TO', 'RD', 'RWD', 'TE'};
    xticklabels(xticklabels_universe)
    title(num2str(isesh))
    
end



%% correlation between sessions using mutually active cells
% all cells is identical because correlations with a vector of zeros are
% always NaN

% preallocate
common_cell_corrs = cell(size(full_session_mtx_3d,3), size(full_session_mtx_3d,3));
common_cell_corrs_cellids = cell(size(common_cell_corrs));
common_cell_corrs_means = nan(size(full_session_mtx_3d,3), size(full_session_mtx_3d,3));
common_cell_matrices = cell(size(full_session_mtx_3d,3)-1,2);

% iterate through session comparisons
event_bins = cumsum([0 time_series_event_spacing]).*100;
for isesh01 = 1:size(full_session_mtx_3d,3)
    for isesh02 = 1:size(full_session_mtx_3d,3)
        
        % active session comparison
        session_mtx01 = full_session_mtx_3d(:,:, isesh01);
        session_mtx02 = full_session_mtx_3d(:,:, isesh02);
        
        % mutual cell idx
        mc_idx = sum(session_mtx01,2)>0 & sum(session_mtx02,2)>0;
        
        % mutual cells
        mc_session01 = session_mtx01(mc_idx,:);
        mc_session02 = session_mtx02(mc_idx,:);
        
        % load common cell matrices
        if isesh02 == isesh01+1
            common_cell_matrices{isesh01, 1} = mc_session01;
            common_cell_matrices{isesh01, 2} = mc_session02;
        end
        
        % only from NP to end of RD
        mc_session01 = mc_session01(:,2:event_bins(5));
        mc_session02 = mc_session02(:,2:event_bins(5));
        
        % iterate through cells
        cell_corrs = nan(size(mc_session01,1), 1);
        for icell = 1:size(mc_session01, 1)
            cell_corrs(icell) = corr(mc_session01(icell,:)', mc_session02(icell,:)');
        end
        
        % load overall
        common_cell_corrs{isesh01, isesh02} = cell_corrs;
        common_cell_corrs_cellids{isesh01, isesh02} = find(mc_idx==1);
        common_cell_corrs_means(isesh01, isesh02) = nanmean(cell_corrs);
        
    end
    
end

% plot means
figure;
ampm_pcolor(common_cell_corrs_means)
axis square
caxis([-1 1])
colorbar
xlabel('Probe'); ylabel('Probe')
title(remove_underscore('cell_corrs_mutual_means'))

%% errorbar line plot of probe over probe similarity

% just probe over probe
cell_corrs_mutual_temp = common_cell_corrs(2:end,1:end-1); 
cell_corrs_mutual_temp = cell_corrs_mutual_temp(logical(eye(size(cell_corrs_mutual_temp,1))));

% transform
for icell = 1:length(cell_corrs_mutual_temp)
    cell_corrs_mutual_temp{icell} = atanh(cell_corrs_mutual_temp{icell});
end

% plot
figure; errorbar_plot(cell_corrs_mutual_temp);

% aesthetics
hold on; plot(xlim, [0 0], 'k--')
ylim([-1 1])
maxr = .99; tic_vect = [-maxr -.9 -.8 -.6 -.3 0 .3 .6 .8 .9 maxr];
ylim(atanh([-maxr maxr])); 
ylim_hold = ylim;
ylim([ylim_hold(1) ylim_hold(2)+(ylim_hold(2)*.25)])
yticks(atanh(tic_vect)); yticklabels(tic_vect)



%% sorting heat maps of adjacent sessions

zoom_portion = 20:430; %end np to end fixed delay

% iterate through comparisons
for icomp = 1:size(common_cell_matrices,1)
    
    
    if size(common_cell_matrices{icomp,1},1)<1
        continue
    end
    
    % sort first session
    [~, sort_idx] = sort_rows_by_peak(common_cell_matrices{icomp,1}(:,zoom_portion));
    
    % plot first session full
    figure; 
    subplot(2,2,1)
    hold on
    imagesc(common_cell_matrices{icomp,1}(sort_idx,:))
    set(gca,'TickLength',[0, 0]); box off;
    title([num2str(icomp) ' full'])
    ylabel('Neuron')
    xlabel('Time')
    
        % red line
        event_frame = cumsum([0 time_series_event_spacing]).*100;
        xticks_hold = [];
        for i = event_frame 
           plot(i.*[1 1], ylim, 'r-', 'linewidth', 1.0) 
           xticks_hold = [xticks_hold i];
        end
        
        % aesthetics
        set(gca, 'YDir','reverse')
        ylim([0+0.5 size(common_cell_matrices{icomp,1}(sort_idx,:),1)+0.5])
        xlim([0.5 size(common_cell_matrices{icomp,1}(sort_idx,:),2)+0.5])
        xticks(xticks_hold)
        xticklabels_universe = {'NP1', 'NP2', 'HE1', 'TO', 'RD', 'RWD', 'TE'};
        xticklabels(xticklabels_universe)
    
    % plot first session zoom
    subplot(2,2,3)
    hold on
    imagesc(common_cell_matrices{icomp,1}(sort_idx, 1:zoom_portion(end)))
    set(gca,'TickLength',[0, 0]); box off;
    title([num2str(icomp) ' zoom'])
    ylabel('Neuron')
    xlabel('Time')
    
        % red line
        event_frame = cumsum([0 time_series_event_spacing]).*100;
        xticks_hold = [];
        for i = event_frame 
           plot(i.*[1 1], ylim, 'r-', 'linewidth', 1.0) 
           xticks_hold = [xticks_hold i];
        end
        
        % aesthetics
        set(gca, 'YDir','reverse')
        ylim([0+0.5 size(common_cell_matrices{icomp,1}(sort_idx, 1:zoom_portion(end)),1)+0.5])
        xlim([0.5 size(common_cell_matrices{icomp,1}(sort_idx, 1:zoom_portion(end)),2)+0.5])
        xticks(xticks_hold)
        xticklabels_universe = {'NP1', 'NP2', 'HE1', 'TO', 'RD', 'RWD', 'TE'};
        xticklabels(xticklabels_universe)
        
    
    
    %sort and plot second session full
    subplot(2,2,2)
    hold on
    imagesc(common_cell_matrices{icomp,2}(sort_idx,:))
    set(gca,'TickLength',[0, 0]); box off;
    title([num2str(icomp+1) ' full'])
    ylabel('Neuron')
    xlabel('Time')
    
        % red line
        event_frame = cumsum([0 time_series_event_spacing]).*100;
        xticks_hold = [];
        for i = event_frame 
           plot(i.*[1 1], ylim, 'r-', 'linewidth', 1.0) 
           xticks_hold = [xticks_hold i];
        end
        
        % aesthetics
        set(gca, 'YDir','reverse')
        ylim([0+0.5 size(common_cell_matrices{icomp,2}(sort_idx,:),1)+0.5])
        xlim([0.5 size(common_cell_matrices{icomp,2}(sort_idx,:),2)+0.5])
        xticks(xticks_hold)
        xticklabels_universe = {'NP1', 'NP2', 'HE1', 'TO', 'RD', 'RWD', 'TE'};
        xticklabels(xticklabels_universe)
    
    %sort and plot second session zoom
    subplot(2,2,4)
    hold on
    imagesc(common_cell_matrices{icomp,2}(sort_idx, 1:zoom_portion(end)))
    set(gca,'TickLength',[0, 0]); box off;
    title([num2str(icomp+1) ' zoom'])
    ylabel('Neuron')
    xlabel('Time')
    
        % red line
        event_frame = cumsum([0 time_series_event_spacing]).*100;
        xticks_hold = [];
        for i = event_frame 
           plot(i.*[1 1], ylim, 'r-', 'linewidth', 1.0) 
           xticks_hold = [xticks_hold i];
        end
        
        % aesthetics
        set(gca, 'YDir','reverse')
        ylim([0+0.5 size(common_cell_matrices{icomp,2}(sort_idx, 1:zoom_portion(end)),1)+0.5])
        xlim([0.5 size(common_cell_matrices{icomp,2}(sort_idx, 1:zoom_portion(end)),2)+0.5])
        xticks(xticks_hold)
        xticklabels_universe = {'NP1', 'NP2', 'HE1', 'TO', 'RD', 'RWD', 'TE'};
        xticklabels(xticklabels_universe)

    set(gcf, 'Position', [219    18   894   863])
end


%% probes only
%{
probe_nums = [1:3:19 20]; 
probes_only = cell(8,8); 
for iprobe1 = 1:8; 
    for iprobe2 = 1:8; 
        probes_only{iprobe1, iprobe2} = common_cell_corrs{probe_nums(iprobe1), probe_nums(iprobe2)}; 
    end; 
end
probes_only_adj = cell(1,8); 
for iprobe1 = 1:7; 
    probes_only_adj{iprobe1} = probes_only{iprobe1,iprobe1+1}; 
end
%}




%% average activity of cells that are held across sessions?

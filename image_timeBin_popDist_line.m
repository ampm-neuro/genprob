function [all_distances, mean_diff_distance, session_mtx_cell_means] = image_timeBin_popDist_line(sesh_num, comp_sesh, trial_num, activity_mtx, cell_regist_mtx, session_number_idx, trial_number_idx, time_bin_idx, time_series_event_spacing)
% plots heatmap of trial activity for all cells, and corresponding lines 
% indicating instantaneous population distance to each comparison session 
% mean

% build session comparison index to restrict neurons
session_comp_idx = cell_regist_mtx(:,sesh_num)>0;
for icomp = 1:length(comp_sesh)
    session_comp_idx = session_comp_idx | cell_regist_mtx(:,comp_sesh(icomp))>0;
end
%session_comp_idx = 1:size(activity_mtx,2);
sesh_trial_act = activity_mtx(session_number_idx==sesh_num & trial_number_idx==trial_num, session_comp_idx)'; 
[mtx, sort_idx] = sort_rows_by_peak(sesh_trial_act);


%% plot sorted firing rates on trial

figure
subplot(2, 3, 1:3); hold on
imagesc(mtx(sum(mtx,2)>0, :))
set(gca,'TickLength',[0, 0]); box off;
title(['Session ' num2str(sesh_num) ', trial number ' num2str(trial_num)])

% red line
event_frame = cumsum([0 time_series_event_spacing]).*100;
xticks_hold = [];
for irl = event_frame(2:4)
   plot(irl.*[1 1], ylim, 'r-', 'linewidth', 1.0) 
end
xlim([0.5 size(mtx(sum(mtx,2)>0, :),2)+0.5])
ylim([0.5 size(mtx(sum(mtx,2)>0, :),1)+0.5])
set(gca, 'YDir','reverse')


%% compute session means
session_mtx_cell_means = nan(length(comp_sesh), sum(session_comp_idx));
for icomp = 1:length(comp_sesh)
    
    % constrain which time bins to include?
    %incl_time_bins = unique(time_bin_idx);
    incl_time_bins = 1:event_frame(4);
    
    % average cell activity
    sesh_idx = session_number_idx==comp_sesh(icomp);
    time_idx = ismember(time_bin_idx, incl_time_bins);
    session_mtx_cell_means(icomp,:) = mean(activity_mtx(sesh_idx & time_idx, session_comp_idx), 1)';
    
end



%% compute instantaneous distances

all_distances = nan(size(sesh_trial_act,2),length(comp_sesh));
for icomp = 1:length(comp_sesh)
    current_comp_sesh = comp_sesh(icomp);
    for ibin = 1:size(sesh_trial_act,2)
        inst_bin_pos = [sesh_trial_act(:,ibin)'; session_mtx_cell_means(icomp,:)];
        all_distances(ibin, icomp) = pdist(inst_bin_pos)./sqrt(size(inst_bin_pos,2));
    end
end



%% plot distances

% colors
colors = parula(length(comp_sesh)+1);
colors(1:end-1,:);

if length(comp_sesh)==2
    subplot(8, 3, 13:18)
else
    subplot(2,3, 4:6)
end

hold on
for icomp = 1:size(comp_sesh,2)
    plot(all_distances(:,icomp), 'color', colors(icomp,:))
end

% red line
for irl = event_frame(2:4)
   plot(irl.*[1 1], ylim, 'r-', 'linewidth', 1.0) 
end

set(gca,'TickLength',[0, 0]); box off;
xlim([1 length(all_distances)])
xticks([])

% legend
legend_trick(colors, '-')
legend_str = cell(1,length(comp_sesh));
for icomp = 1:length(comp_sesh)
    legend_str{icomp} = num2str(comp_sesh(icomp));
end
legend(legend_str)
set(gca,'TickLength',[0, 0]); box off;





%% compute and plot differences

% preallocate to avoid output errors if number of comparisons ~=2
mean_diff_distance = [];

%dotLine_rng = incl_time_bins(end)+1:size(all_distances,1);
dotLine_rng = event_frame(3):event_frame(5);

% plot
if length(comp_sesh)==2

    subplot(8, 3, 19:24)
    hold on

    plot(all_distances(:,1)-all_distances(:,2), 'k')
    plot(xlim, [1 1].*mean(all_distances(dotLine_rng,1)-all_distances(dotLine_rng,2)), 'k--')
    xlim([1 length(all_distances)])
    set(gca,'TickLength',[0, 0]); box off;
    
    mean_diff_distance = mean(all_distances(dotLine_rng,1)-all_distances(dotLine_rng,2));
    
    % red line
    for irl = event_frame(2:4)
       plot(irl.*[1 1], ylim, 'r-', 'linewidth', 1.0) 
    end

elseif length(comp_sesh)==1

    mean_diff_distance = mean(all_distances(dotLine_rng));

end




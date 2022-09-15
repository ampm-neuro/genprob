function ALL_GenLearn(day_bin_size, accordian_bin_size)
%runs GenLearn for multiple learning stages

%figures
h1 = figure; hold on;
set(gca, 'XScale', 'log')
set(gca,'TickLength',[0, 0]); box off;
h2 = figure; hold on;
set(gca,'TickLength',[0, 0]); box off;

%max number of sessions
max_sesh = 0;
folderpath = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\Gen\';
file_list_subjects = dir(folderpath);
file_list_subjects(1:2) = [];
for isubject = 1:length(file_list_subjects)
    file_list_sessions = dir([folderpath '\' file_list_subjects(isubject).name]);
    max_sesh = max([max_sesh length(file_list_sessions(3:end))]);
end
max_sesh_original = max_sesh;

%divide into rows for each day_bin
if rem(max_sesh, day_bin_size)
    max_sesh = max_sesh + (day_bin_size - rem(max_sesh, day_bin_size));
end
cols = max_sesh/day_bin_size;
day_bins = reshape(1:max_sesh, day_bin_size, cols)';

%colors and legend hack
colors = distinguishable_colors(size(day_bins,1));
figure(h1); for i = 1:size(colors,1); plot(15000,0,'markersize', realmin, 'color', colors(i,:)); end

% preallocate
day_bin_waits = cell(1, size(day_bins,1));
day_bin_fielddiff_zs = cell(1);
z_means = nan(size(day_bins,1),1);
z_ses = nan(size(day_bins,1),1);
max_col = 0;
unique_subjs = [];

for day_bin = 1:size(day_bins,1)
    zdiffs_hold = [];
    
    for day = day_bins(day_bin,:)

        % catch uncollected day
        if day > max_sesh_original
            continue
        end
        
        % compute waits for day
        [wait_durations_day, unq_frq, p_dist, zdiffs] = GenLearn(day);

        %{
        if day == 9
           for i = 1:length(wait_durations_day)
               wait_durations_day{i}
                wait_durations_day{i} = wait_durations_day{i}(1);
                wait_durations_day{i}
           end
        end
        %}
        
        unique_subjs = unique([unique_subjs; zdiffs(:,1)]);
        
        %accomodate missing subject days
        [zdiffs_hold, mc] = match_subjects(zdiffs_hold, zdiffs);
        max_col = max([max_col mc]);
        
        % concatenate day waits into day bin
        if day == min(day_bins(day_bin,:))
            day_bin_waits{day_bin} = wait_durations_day;
        else
            for ifreq = 1:length(day_bin_waits{day_bin})               
                day_bin_waits{day_bin}{ifreq} = [day_bin_waits{day_bin}{ifreq}; wait_durations_day{ifreq}];
            end
        end

        day_bin_fielddiff_zs{day_bin} = nanmean(zdiffs_hold,3);
        z_means(day_bin) = nanmean(day_bin_fielddiff_zs{day_bin}(:,2),1);
        z_ses(day_bin) = nanstd(day_bin_fielddiff_zs{day_bin}(:,2),1)./sqrt(sum(~isnan(day_bin_fielddiff_zs{day_bin}(:,2)),1));
        
    end
    
    
    % bin tones
    %
    
    if rem(length(unq_frq),accordian_bin_size)~=0
        error('bin size must evenly divide number of tones'); 
    end

    se_wait_durations = nan(1,length(day_bin_waits{day_bin})/accordian_bin_size);
    mean_wait_durations = nan(size(se_wait_durations));
    count = 0;
    for ise = reshape(1:length(unq_frq), accordian_bin_size, length(unq_frq)/accordian_bin_size)
        count = count+1;

        wait_durations_hold = cell2mat(day_bin_waits{day_bin}(ise));

        if length(wait_durations_hold)>=5     
            se_wait_durations(count) = nanstd(wait_durations_hold)./sqrt(sum(~isnan(wait_durations_hold)));
            mean_wait_durations(count) = nanmean(wait_durations_hold);
        else
            se_wait_durations(count) = nan; 
            mean_wait_durations(count) = nan;
        end
    end
    unq_frq_plot = mean(reshape(unq_frq, accordian_bin_size, length(unq_frq)/accordian_bin_size),1);
    p_dist_plot = mean(reshape(p_dist, accordian_bin_size, length(p_dist)/accordian_bin_size),1);

    %plots
    %
    figure(h1)
    %errorbar(unq_frq_plot, mean_wait_durations, se_wait_durations, 'linewidth', 1.5)
    plot(unq_frq_plot, smooth(inpaint_nans(mean_wait_durations), 3), '-', 'linewidth', 1, 'color', colors(day_bin,:));
    plot(unq_frq_plot, smooth(inpaint_nans(mean_wait_durations) + inpaint_nans(se_wait_durations), 3), '-', 'linewidth', 0.5, 'color', colors(day_bin,:));
    plot(unq_frq_plot, smooth(inpaint_nans(mean_wait_durations) - inpaint_nans(se_wait_durations), 3), '-', 'linewidth', 0.5, 'color', colors(day_bin,:));

    
end

%match subjects across zdists
for iday = 1:length(day_bin_fielddiff_zs)
    if size(day_bin_fielddiff_zs{iday},1) < max_col
        nan_hold = nan(max_col, size(day_bin_fielddiff_zs{iday},2), size(day_bin_fielddiff_zs{iday},3));
        nan_hold(day_bin_fielddiff_zs{iday}(:,1),:, :) = day_bin_fielddiff_zs{iday};
        day_bin_fielddiff_zs{iday} = nan_hold;
    end
    day_bin_fielddiff_zs{iday} =  day_bin_fielddiff_zs{iday}(:,2,:);
end

%plot aesthetics
figure(h1)
xlim_hold = xlim;
plot(xlim, [0 0], 'k--', 'linewidth', 4)
plot([min(unq_frq_plot(p_dist_plot>0.5)) max(unq_frq_plot(p_dist_plot>0.5))],[0 0], 'r-', 'linewidth', 4)
ylabel('Wait times')
xlabel('Tone Frequency')
legend_input = cell(1);
for istage = 1:size(day_bins,1)
    legend_input{istage} = ['Stage ' num2str(istage)];
end
legend(legend_input, 'location','northeastoutside')

figure(h2)
%plot(1:day_bin, cell2mat(day_bin_fielddiff_zs)', '-o', 'color', 0.8.*[1 1 1])
plot(1:day_bin, cell2mat(day_bin_fielddiff_zs)', '-o')
legend
errorbar(z_means, z_ses, 'k-', 'linewidth', 1.5)
xlim([.5 day_bin+.5]); xticks(1:day_bin); xlabel('Training stage')
ylim([-5 5])
hold on; plot(xlim, [0 0], 'k--')
for istage = 1:day_bin
    if ttest(day_bin_fielddiff_zs{istage})
        text(istage,4,'*','FontSize',30,'HorizontalAlignment','center')
    end
end
end


function [mtx_out, max_col] = match_subjects(mtx1, mtx2)
%takes two matrices and 3d concanonates them by matching numbers in
%first column. Unique items become nan rows in other mtx.
    
if isempty(mtx1)
    mtx_out = mtx2;
    max_col = max(mtx2(:,1));
    return
elseif isempty(mtx2)
    mtx_out = mtx1;
    max_col = max(mtx1(:,1));
    return
end

max_col = max([mtx1(:,1);mtx2(:,1)]);

mtx_out = nan(max_col,2,size(mtx1,3)+size(mtx2,3));

mtx_out(mtx1(:,1),2,1:size(mtx1,3)) = repmat(mtx1(:,2), 1, 1, size(mtx1,3));
mtx_out(mtx2(:,1),2,(size(mtx1,3)+1): (size(mtx1,3)+size(mtx2,3))) = repmat(mtx2(:,2), 1, 1, size(mtx2,3));

mtx_out(:,1,:) = repmat(1:max_col,1,1,size(mtx_out,3));

end



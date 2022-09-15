function [all_in, all_out] = ALL_betweensesh_learn(stage_num)

%should be input
folderpath = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\Gen_richards';

%preallocate
day_hold_in_means = cell(1,100);
day_hold_out_means = cell(1,100);
day_hold_in_all = cell(1,100);
day_hold_out_all = cell(1,100);
zdists_all = [];%cell(1,100);
legend_input = cell(1);
all_in = [];
all_out = [];

%figure
h1 = figure; hold on
h2 = figure; hold on
xmax = 0;

%subject folders
file_list_subjects = dir(folderpath);
file_list_subjects(1:2) = [];

%cage manipulations
incl_subj = [1:4];
m_per_cage = 4;
all_subjs = 1:length(file_list_subjects);
all_subs = reshape(all_subjs, length(all_subjs)/m_per_cage, m_per_cage);
all_subs = unique(all_subs(incl_subj, :));
if size(all_subs,1)>size(all_subs,2)
    all_subs = all_subs';
end
%all_subs=all_subs(1:4:end)
%all_subs = [1:4];
%all_subs = [5:8];
%all_subs = [9:12];
%all_subs = [13:16];

%colors
colors = distinguishable_colors(length(all_subs));

subj_count = 0;
for isubject = all_subs%1:length(file_list_subjects)
    current_subj = file_list_subjects(isubject).name;
    subj_count = subj_count +1;
    
    %session files
    file_list_sessions = dir([folderpath '\' current_subj]);
    file_list_sessions(1:2) = [];  
    
    %preallocate
    subj_hold_means = [];
    subj_hold_rich_all = cell(1,2);
    subj_hold_poor_all = cell(1,2);
    included_session_ct = 0;

    
    for isession = 1 : length(file_list_sessions)
        
        current_sesh = file_list_sessions(isession).name;
        
        %skip if session file does not exist
        if length(file_list_sessions)<isession
            continue
        elseif ~contains(file_list_sessions(isession).name, 'd')
        end
        
        %skip sessions from other learning stages
        stage_num = num2str(stage_num);
        if length(stage_num)<2
           stage_num = ['0' stage_num]; 
        end
        if ~contains(current_sesh, [ 'train' stage_num])
           continue 
        end

        %load data file
        current_sesh = file_list_sessions(isession).name;
        load([folderpath '\' current_subj '\' current_sesh])
        
        
        
        %compute wait durations
        [mean_wait_durations, wait_durations, ~, unq_frq, p_dist] = wait_times(trl_mtx,medass_cell,0);
        
        
        %trial minimum
        rich_min = 2;
        if length(wait_durations{1})<rich_min
            continue
        end
        included_session_ct = included_session_ct + 1;
        
        
        
        %load output
        subj_hold_means = [subj_hold_means; [nanmean(mean_wait_durations(p_dist>.5)) nanmean(mean_wait_durations(p_dist<.5))]];
        subj_hold_rich_all{included_session_ct} = cell2mat(wait_durations(p_dist>0.50));
        subj_hold_poor_all{included_session_ct} = cell2mat(wait_durations(p_dist<0.50));
        day_hold_in_means{included_session_ct} = [day_hold_in_means{included_session_ct};nanmean(mean_wait_durations(p_dist>.5))];
        day_hold_out_means{included_session_ct} = [day_hold_out_means{included_session_ct};nanmean(mean_wait_durations(p_dist<.5))];
        day_hold_in_all{included_session_ct} = [day_hold_in_all{included_session_ct};cell2mat(wait_durations(p_dist>0.50))];
        day_hold_out_all{included_session_ct} = [day_hold_out_all{included_session_ct};cell2mat(wait_durations(p_dist<0.50))];
        

    end

    %no included sessions
    if isempty(subj_hold_means)
        continue
    end
    
    xmax = max([xmax included_session_ct]);
    
    figure(h1)
    %infield
    plot(1:size(subj_hold_means,1), subj_hold_means(:,1), '--o', 'color', colors(subj_count,:), 'linewidth', 1)
    %outfield
    plot(1:size(subj_hold_means,1), subj_hold_means(:,2), '-o', 'color', colors(subj_count,:), 'linewidth', 1)
    
    figure(h2)
    %infield-outfield
    %plot(1:size(subj_hold,1), subj_hold(:,1)-subj_hold(:,2), '-o', 'color', [0 0 0], 'linewidth', 1)
    means_z=nan(2,length(subj_hold_rich_all));stds_z=nan(size(means_z));zdists=nan(1,length(subj_hold_rich_all));
    for iz = 1:length(subj_hold_rich_all)
        means_z(:,iz) = [mean(subj_hold_rich_all{iz}); mean(subj_hold_poor_all{iz})];
        stds_z(:,iz) = [std(subj_hold_rich_all{iz}); std(subj_hold_poor_all{iz})];
        zdists(iz) = (means_z(1,iz)-means_z(2,iz))/mean(stds_z(:,iz));
    end

    
    %accomodate inequal number of training sessions
    if size(zdists,2) < size(zdists_all,2)
        zdists = [zdists nan(size(zdists,1),size(zdists_all,2)-size(zdists,2))];
    elseif size(zdists,2) > size(zdists_all,2)
        zdists_all = [zdists_all nan(size(zdists_all,1), size(zdists,2)-size(zdists_all,2))];
    elseif ~isempty(zdists_all) && size(zdists,2) < size(zdists_all,2)
        zdists_all = [zdists_all nan(size(zdists,2)-size(zdists_all,2),size(zdists_all,1))];
    end
    
    %plot and plot prep
    zdists_all = [zdists_all; zdists];
    plot(1:size(zdists,2), zdists, '-o', 'color', colors(subj_count,:), 'linewidth', 1)
    legend_input =[legend_input {current_subj}];
    
end
%plot error bar
errorbar(1:size(zdists_all,2), nanmean(zdists_all), nanstd(zdists_all)./sqrt(sum(~isnan(zdists_all))), 'k', 'linewidth', 2)

legend(legend_input(2:end), 'location', 'NorthEastOutside')

%add pvals
for ip = 1:size(zdists_all,2)

    [~, pval_z] =  ttest( zdists_all(:, ip) , 0);
    
   if ~isempty(pval_z)
      hold on; text( ip-.15, 4.0, num2str(pval_z))
   end
end

%plot means
figure(h1)
means_in = nan(1, length(day_hold_in_means));
means_out = nan(1, length(day_hold_out_means));
ses_in = nan(1, length(day_hold_in_means));
ses_out = nan(1, length(day_hold_out_means));
for ix = 1:length(day_hold_in_means)
    means_in(ix) = mean(day_hold_in_means{ix});
    means_out(ix) = mean(day_hold_out_means{ix});
    ses_in(ix) = std(day_hold_in_means{ix})/sqrt(length(day_hold_in_means{ix}));
    ses_out(ix) = std(day_hold_out_means{ix})/sqrt(length(day_hold_out_means{ix}));
end
errorbar(means_in, ses_in, '--', 'color', [0 0 0], 'linewidth', 4)
errorbar(means_out, ses_out, '-', 'color', [0 0 0], 'linewidth', 4)


%all trials figure
celle = cell(2,1);
f1 = figure; hold on
f2 = figure; hold on
for i = 1:xmax
    
    figure(f1)
    subplot(1,xmax, i);
    celle{2} = day_hold_in_all{i}; celle{1} = day_hold_out_all{i};
    errorbar_plot(celle, 0)
    ylim([-2 38])
    xticks([1 2])
    xticklabels({'Poor tone(s)','Rich tone(s)'})
    yticks(-2:5:38)
    yticklabels(yticks+2)
    
    figure(f2)
    subplot(xmax, 1, i);
    hold on
    histogram(day_hold_out_all{i}, 0:1:40, 'normalization', 'probability')
    histogram(day_hold_in_all{i}, 0:1:40, 'normalization', 'probability')
    legend({'Poor tone', 'Rich tone'})
    set(gca,'TickLength',[0, 0]); box off;
    ylim([0 .3])
    xlim([0 40])
    
    all_in = [all_in; day_hold_in_all{i}];
    all_out = [all_out; day_hold_out_all{i}];
end



%beautification
figure(h1);
xlim([0.5 xmax+0.5])
xticks([1:1:xmax])
ylim([0 inf])
set(gca,'TickLength',[0, 0]);
legend({'Rich', 'Poor'}, 'location', 'NorthEastOutside')
xlabel('Training day')
ylabel('Mean wait time (s)')
%xlim([4500 42000])
%xticks([5000 8500 14000 23000 35000])

figure(h2);
xlim([0.5 xmax+0.5])
xticks([1:1:xmax])
ylim([-2.5 4.25])
plot(xlim,[0 0], 'k--')
set(gca,'TickLength',[0, 0]);
xlabel('Training day')
ylabel('Rich - Poor (pooled var)')
%xlim([4500 42000])
%xticks([5000 8500 14000 23000 35000])


function [all_stage_learn_mtx, all_stage_learn_cell, days_to_crit, all_TwoDay_learn_mtx] = ALL_stage_learn(datafolder, stages, num_stage_samples)
%plots learning curves from each training stage for all animals. input
%dictates how to divide up the sessions. 
%

% could be input
folderpath = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\';
folderpath = [folderpath datafolder];

% file list
file_list_subjects = dir(folderpath);
file_list_subjects(1:2) = [];

%subjects
num_subj = length(file_list_subjects);

% colors
colors = distinguishable_colors(num_subj);

% preallocate combined output
all_stage_learn_cell = cell(1,length(stages));
all_stage_learn_mtx = nan(num_subj, num_stage_samples, length(stages));
all_TwoDay_learn_mtx = nan(num_subj, num_stage_samples, length(stages));
num_rich_sessions = nan(num_subj, num_stage_samples, length(stages));

%% iterate through stages
stage_count = 0;
for istage = stages
    stage_count = stage_count + 1;
    
    % compute learning curves for all subjects
    [subject_cell, rwd_bias_subj_cell] = ALL_stage_LearnCurves(folderpath, istage);
        
    % remove sessions with too few probe trials
    min_rich_probe_trials = 3; % minimum samples
    for isubj = 1:num_subj

       delete_idx = [];
       for isesh = 1:size(subject_cell{isubj}, 2)

            %check number of rich probes
            if length(subject_cell{isubj}{1, isesh})<min_rich_probe_trials
                delete_idx = [delete_idx isesh];
            end
       end
       
       %delete sessions with too few samples
       subject_cell{isubj}(:,delete_idx) = [];
       
    end


    % preallocate stage summaries
    stage_learn_cell = cell(num_subj, num_stage_samples);
    stage_learn_mtx = nan(num_subj, num_stage_samples);
    stage_learn_mtx_2day = nan(num_subj, 2);
    stage_rwd_bias_mtx = nan(num_subj, num_stage_samples);

    % iterate through subjects
    for isubj = 1:num_subj

        num_sessions = size(subject_cell{isubj},2);
        if num_sessions == 0
            continue
        end

            target_sessions = round(linspace(1, num_sessions, num_stage_samples));
            stage_sample_count = 0;
            for stage_sample = 1:length(target_sessions)
                stage_sample_count = stage_sample_count + 1;

                % rich and poor tone wait times
                rich = subject_cell{isubj}{1, target_sessions(stage_sample)};
                poor = subject_cell{isubj}{2, target_sessions(stage_sample)};

                % difference between means in terms of std
                dprime_hold = computeCohen_d(rich, poor);
                
                % load rich and poor
                stage_learn_cell{isubj, stage_sample_count} = [{rich} {poor}];
                stage_learn_mtx(isubj, stage_sample_count) = dprime_hold;
                
                %load recency bias
                stage_rwd_bias_mtx(isubj, stage_sample_count) = rwd_bias_subj_cell{isubj}{1, target_sessions(stage_sample)};
            end
            
            
            % second day
            if num_sessions >= 2
                dprime_hold = computeCohen_d(subject_cell{isubj}{1, 2}, subject_cell{isubj}{2, 2});
                stage_learn_mtx_2day(isubj, 2) = dprime_hold;
            end
    end
    
    % add first day to 2day mtx
    stage_learn_mtx_2day(:,1) = stage_learn_mtx(:,1);
    
    %load stage summaries
    all_stage_learn_cell{stage_count} = stage_learn_cell;
    all_stage_learn_mtx(:, :, stage_count) = stage_learn_mtx;
    all_TwoDay_learn_mtx(:, :, stage_count) = stage_learn_mtx_2day;
    
end

for i = 2:6
    %[a b c d] = ttest2(all_stage_learn_mtx(:,1,1), all_stage_learn_mtx(:,1,i))
end

%% stage curve figure
figure; hold on
stage_firsts_zs = nan(num_subj, length(stages));
stage_lasts_zs = nan(num_subj, length(stages));
stage_firsts_rwd_bias = nan(num_subj, length(stages));
for istage = 1:size(all_stage_learn_mtx,3)
    
    %subplot(length(stages), 1, istage)
    rng = 1/3;
    xpos = linspace(istage-rng,istage+rng, num_stage_samples);
    
    for isubj = 1:num_subj
       plot(xpos, all_stage_learn_mtx(isubj,:,istage), 'o-', 'color', colors(isubj,:))

       %capture first and last day performances
       stage_firsts_zs(isubj,istage) = all_stage_learn_mtx(isubj,1,istage);
       stage_lasts_zs(isubj,istage) = all_stage_learn_mtx(isubj,end,istage);
       
    end
    
    errorbar(xpos, nanmean(all_stage_learn_mtx(:,:,istage)), ...
        nanstd(all_stage_learn_mtx(:,:,istage))...
        ./sqrt(sum(~isnan(all_stage_learn_mtx(:,:,istage)))), 'k-', 'linewidth', 3.5)

end
% one legend
legend(file_list_subjects.name, 'location', 'northeastoutside')
set(gca,'TickLength',[0, 0]); box off;
ylim([-2.5 5])
xlim([min(stages)-0.75 max(stages)+0.75])
plot(xlim, [0 0], 'k--')
xlabel('Training stage')
ylabel('Preference for rich tone')


%% plot first day performance (z)
figure; hold on;
for isubj = 1:num_subj
    plot(1:size(stage_firsts_zs,2), stage_firsts_zs(isubj,:), 'o-', 'color', colors(isubj,:))
end
errorbar(nanmean(stage_firsts_zs), ...
        nanstd(stage_firsts_zs)./sqrt(sum(~isnan(stage_firsts_zs))), 'k-', 'linewidth', 3.5);
title('First day zs')
legend(file_list_subjects.name, 'location', 'northeastoutside')
set(gca,'TickLength',[0, 0]); box off;
xlim([min(stages)-0.5 max(stages)+0.5])
hold on; plot(xlim, [1 1].*0, 'k--')

    % bar
    figure; hold on
    bar([nanmean(stage_firsts_zs(:,1)), nanmean(stage_firsts_zs(:,end))])
    errorbar_plot([{stage_firsts_zs(:,1)} {stage_firsts_zs(:,end)}], 1)
    ylabel('dprime; first day perform')



%% plot last day performance minus first day performance
figure; hold on;
for isubj = 1:num_subj
    plot(1:size(stage_lasts_zs - stage_firsts_zs,2), stage_lasts_zs(isubj,:) - stage_firsts_zs(isubj,:), 'o-', 'color', colors(isubj,:))
end
errorbar(nanmean(stage_lasts_zs - stage_firsts_zs), ...
        nanstd(stage_lasts_zs - stage_firsts_zs)./sqrt(sum(~isnan(stage_lasts_zs - stage_firsts_zs))), 'k-', 'linewidth', 3.5);
title('Last day zs - first day zs')
legend(file_list_subjects.name, 'location', 'northeastoutside')
set(gca,'TickLength',[0, 0]); box off;
xlim([min(stages)-0.5 max(stages)+0.5])
hold on; plot(xlim, [1 1].*0, 'k--')

    % bar
    figure; hold on
    bar([nanmean(stage_lasts_zs(:,1) - stage_firsts_zs(:,1)), nanmean(stage_lasts_zs(:,end) - stage_firsts_zs(:,end))])
    errorbar_plot([{stage_lasts_zs(:,1) - stage_firsts_zs(:,1)} {stage_lasts_zs(:,end) - stage_firsts_zs(:,end)}], 1)
    ylabel('dprime; amount learning')


%% days to criterion line graph and histograms

% compute all days to criterion
days_to_crit = nan(num_subj, length(stages));
for istage = 1:length(stages)
    subject_cell = ALL_stage_LearnCurves(folderpath, istage);
    for isubj = 1:num_subj
        days_to_crit(isubj,istage) = 0;
        for isesh=1:size(subject_cell{isubj},2)
               days_to_crit(isubj,istage) = days_to_crit(isubj,istage) +1;
               
               % break if subject_cell{isubj} meets crit requirements
               [~,pval,~,stats] = ttest2(subject_cell{isubj}{1,isesh},subject_cell{isubj}{2,isesh});
               if pval<0.01 && stats.tstat>0
                   break
               end
        end
    end
end

days_to_crit(days_to_crit==0)=nan;
days_to_crit(days_to_crit>10) = 10;

    % bar
    figure; hold on
    bar([nanmean(days_to_crit(:,1)), nanmean(days_to_crit(:,end))])
    errorbar_plot([{days_to_crit(:,1)} {days_to_crit(:,end)}], 1)
    ylabel('days to crit')
    ylim([0 10.75])

for i = 2:6
    %[a b c d] = ttest2(days_to_crit(:,1), days_to_crit(:,i))
end
    
figure; hold on

%% plot subject values
xpos = 1:size(days_to_crit,2);
for isubj = 1:num_subj
    jx = jitter_xpos(xpos, days_to_crit(isubj,:));
    plot(jx, days_to_crit(isubj,:), '-o', 'color', colors(isubj,:))
end
errorbar(nanmean(days_to_crit), nanstd(days_to_crit)./sqrt(sum(~isnan(days_to_crit))), 'k-', 'linewidth', 3.5);
set(gca,'TickLength',[0, 0]); box off;
ylim([0 max(days_to_crit(~isnan(days_to_crit)))+1])
yticks(0:max(days_to_crit(~isnan(days_to_crit)))+1)
xlim([min(stages)-0.5 max(stages)+0.5])
xticks(min(stages):max(stages))
xlabel('Training stage')
ylabel('Days to crit')
legend(file_list_subjects.name, 'location', 'northeastoutside')




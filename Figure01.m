green_color = [21 154 72]./255; blue_color =  [43 73 159]./255;
light_green_color = [162 201 129]./255; light_blue_color = [102 164 213]./255;
figure_path = 'C:\Users\ampm1\Documents\manuscripts\gentrain\manuscript_figures\fig01\';

%% Discrimination performance
[firstLast_dprime_mtx_predict, ~, days_to_crit_predict] = ALL_stage_learn(['two_tone\' 'train_mevar'], 1:6, 2);
[firstLast_dprime_mtx_unpredict, ~, days_to_crit_unpredict] = ALL_stage_learn(['two_tone\' 'train_hivar'], 1:6, 2);
for iplot = 1:14; close; end % close legacy plots


%% discrimination line plots

    % predict
    figure; hold on
    for istage = 1:size(firstLast_dprime_mtx_predict,3)
        errorbar_plot([{firstLast_dprime_mtx_predict(:,1,istage)} {firstLast_dprime_mtx_predict(:,end,istage)}], 1, [-0.3 +0.3]+istage, light_green_color, green_color)
    end
    xlabel('Training problem (first and last day)')
    ylabel('Preference for rich tone (d`)')
    axis([0.25 6.75 -3 10]); axis square
    plot(xlim, [0 0], 'k--')
    title('Predictable learning curve')
    %var_name = 'LearningCurve_predictable'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
    
        % stats: first days compared to day1
        % rm anova
        rmtbl_firstdays_to_d1_predict = simple_mixed_anova(squeeze(firstLast_dprime_mtx_predict(:,1,:)));
        
        % pairwise
        pvals_firstday_predict_d1 = nan(1,5); stats_cell_firstday_predict_d1 = cell(1,5);
        for itest = 1:5
            [~, pvals_firstday_predict_d1(itest), ~, stats_cell_firstday_predict_d1{itest}] = ttest(firstLast_dprime_mtx_predict(:,1,1), firstLast_dprime_mtx_predict(:,1,itest+1));
        end
        
        % stats: compared to 0 (not super informative)
        pvals_firstday_predict_0 = nan(1,5); stats_cell_firstday_predict_0 = cell(1,6);
        for itest = 1:6
            [~, pvals_firstday_predict_0(itest), ~, stats_cell_firstday_predict_0{itest}] = ttest(firstLast_dprime_mtx_predict(:,1,itest));
        end
    
    % unpredict
    figure; hold on
    for istage = 1:size(firstLast_dprime_mtx_unpredict,3)
        errorbar_plot([{firstLast_dprime_mtx_unpredict(:,1,istage)} {firstLast_dprime_mtx_unpredict(:,end,istage)}], 1, [-0.3 +0.3]+istage, light_blue_color, blue_color)
    end
    xlabel('Training problem (first and last day)')
    ylabel('Preference for rich tone (d`)')
    axis([0.25 6.75 -3 10]); axis square
    plot(xlim, [0 0], 'k--')
	title('Predictable learning curve')
    %var_name = 'LearningCurve_unpredictable'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
    
        % stats: first days compared to day1
        % repeated measures anova
        rmtbl_firstdays_to_d1_unpredict = simple_mixed_anova(squeeze(firstLast_dprime_mtx_unpredict(:,1,:)));
        
        % pairwise ttests
        pvals_firstday_unpredict_d1 = nan(1,5); stats_cell_firstday_unpredict_d1 = cell(1,5);
        for itest = 1:5
            [~, pvals_firstday_unpredict_d1(itest), ~, stats_cell_firstday_unpredict_d1{itest}] = ttest(firstLast_dprime_mtx_unpredict(:,1,1), firstLast_dprime_mtx_unpredict(:,1,itest+1));
        end
        
        % stats: compared to 0 (not super informative)
        pvals_firstday_unpredict_0 = nan(1,5); stats_cell_firstday_unpredict_0 = cell(1,6);
        for itest = 1:6
            [~, pvals_firstday_unpredict_0(itest), ~, stats_cell_firstday_unpredict_0{itest}] = ttest(firstLast_dprime_mtx_unpredict(:,1,itest));
        end

        
        
%% days to criterion line plots
        
    % predict
    figure; hold on
    d2c_cell_predict = cell(1,size(days_to_crit_predict,2));
    for iprob = 1:size(days_to_crit_predict,2)
        d2c_cell_predict{iprob} = days_to_crit_predict(:,iprob);
    end
    errorbar_plot(d2c_cell_predict, 1, 1:6, light_green_color, green_color)
    xlabel('Training problem')
    ylabel('Days to criterion')
    axis([0.75 6.25 .5 10.5]); axis square; xticks(1:6)
    plot(xlim, [0 0], 'k--')
    title('Predictable learning rates')
    %var_name = 'LearningRates_predictable'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
    
        % stats: first days compared to day1
        % repeated measures anova
        rmtbl = simple_mixed_anova(days_to_crit_predict);
        ranova_d2c_stats_predict = [rmtbl.DF(2) rmtbl.DF(3) rmtbl.F(3) rmtbl.pValue(3)];
        
        % pairwise ttests
        pvals_d2c_predict = nan(1,5); stats_cell_d2c_predict_d1 = cell(1,5);
        for itest = 1:5
            [~, pvals_d2c_predict(itest), ~, stats_cell_d2c_predict_d1{itest}] = ttest(days_to_crit_predict(:,1), days_to_crit_predict(:,itest+1));
        end
    
    % predict
    figure; hold on
    d2c_cell_unpredict = cell(1,size(days_to_crit_unpredict,2));
    for iprob = 1:size(days_to_crit_unpredict,2)
        d2c_cell_unpredict{iprob} = days_to_crit_unpredict(:,iprob);
    end
    errorbar_plot(d2c_cell_unpredict, 1, 1:6, light_blue_color, blue_color)
    xlabel('Training problem')
    ylabel('Days to criterion')
    axis([0.75 6.25 .5 10.5]); axis square; xticks(1:6)
    plot(xlim, [0 0], 'k--')
    title('Unpredictable learning rates')
    %var_name = 'LearningRates_unpredictable'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
    
        % stats: first days compared to day1
        % repeated measures anova
        rmtbl = simple_mixed_anova(days_to_crit_unpredict);
        ranova_d2c_stats_unpredict = [rmtbl.DF(2) rmtbl.DF(3) rmtbl.F(3) rmtbl.pValue(3)];
        
        % pairwise ttests
        pvals_d2c_unpredict = nan(1,5); stats_cell_d2c_unpredict_d1 = cell(1,5);
        for itest = 1:5
            [~, pvals_d2c_unpredict(itest), ~, stats_cell_d2c_unpredict_d1{itest}] = ttest(days_to_crit_unpredict(:,1), days_to_crit_unpredict(:,itest+1));
        end
        
%% direct comparison plots

% first day
figure; hold on
errorbar_plot([{firstLast_dprime_mtx_predict(:,1,1)} {firstLast_dprime_mtx_predict(:,1,end)}], 1, [1 2], light_green_color, green_color)
errorbar_plot([{firstLast_dprime_mtx_unpredict(:,1,1)} {firstLast_dprime_mtx_unpredict(:,1,end)}], 1, [1 2], light_blue_color, blue_color)
axis([0.75 2.25 -2 6]); %axis square
plot(xlim, [0 0], 'k--')
xticks(1:2); xticklabels({'1', '6'}); xlabel('Problem')
ylabel('Preference for rich tone (d`)')
title('First day discrimination')
set(gcf, 'Position', [680   558   337   420])
%var_name = 'FirstDayDiscr_compare'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')

    %  stats: mixed model test for interaction
    % repeated measures anova
    datamtx = squeeze([firstLast_dprime_mtx_predict(:, 1, [1 6]); firstLast_dprime_mtx_unpredict(:, 1, [1 6])]);
    between_subjs = [zeros(40,1); ones(40,1)];
    rmtbl = simple_mixed_anova(datamtx, between_subjs, {'problem'}, {'group'});
    
    % ttest2s for group comparisons
    [~, pval_firstday_comp_problem1, ~, stats_firstday_comp_problem1] = ttest2(firstLast_dprime_mtx_predict(:,1,1), firstLast_dprime_mtx_unpredict(:,1,1));
    [~, pval_firstday_comp_problem6, ~, stats_firstday_comp_problem6] = ttest2(firstLast_dprime_mtx_predict(:,1,6), firstLast_dprime_mtx_unpredict(:,1,6));
    

% days to criterion
figure; hold on
errorbar_plot([{days_to_crit_predict(:,1)} {days_to_crit_predict(:,6)}], 1, [1 2], light_green_color, green_color)
errorbar_plot([{days_to_crit_unpredict(:,1)} {days_to_crit_unpredict(:,6)}], 1, [1 2], light_blue_color, blue_color)
axis([0.75 2.25 .5 10.5]); %axis square
plot(xlim, [0 0], 'k--')
xticks(1:2); xticklabels({'1', '6'}); xlabel('Problem')
ylabel('Days to criterion')
title('Learning rate')
set(gcf, 'Position', [680   558   337   420])
%var_name = 'LearningRate_compare'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')


    % stats: mixed model test for interaction
    % repeated measures anova
    datamtx = [days_to_crit_predict(:, [1 6]); days_to_crit_unpredict(:, [1 6])];
    between_subjs = [zeros(40,1); ones(40,1)];
    rmtbl = simple_mixed_anova(datamtx, between_subjs, {'problem'}, {'group'});
    ranova_d2c_stats_comp = [rmtbl.DF(3) rmtbl.DF(5) rmtbl.F(5) rmtbl.pValue(5)];

    % ttest2s for group comparisons
    [~, pval_d2c_comp_problem1, ~, stats_d2c_comp_problem1] = ttest2(days_to_crit_predict(:,1), days_to_crit_unpredict(:,1));
    [~, pval_d2c_comp_problem6, ~, stats_d2c_comp_problem6] = ttest2(days_to_crit_predict(:,6), days_to_crit_unpredict(:,6));
   

% amount of learning
figure; hold on
first_day_diff_predict = firstLast_dprime_mtx_predict(:,2,1) - firstLast_dprime_mtx_predict(:,1,1);
last_day_diff_predict = firstLast_dprime_mtx_predict(:,2,6) - firstLast_dprime_mtx_predict(:,1,6);
first_day_diff_unpredict = firstLast_dprime_mtx_unpredict(:,2,1) - firstLast_dprime_mtx_unpredict(:,1,1);
last_day_diff_unpredict = firstLast_dprime_mtx_unpredict(:,2,6) - firstLast_dprime_mtx_unpredict(:,1,6);
errorbar_plot([{first_day_diff_predict} {last_day_diff_predict}], 1, [1 2], light_green_color, green_color)
errorbar_plot([{first_day_diff_unpredict} {last_day_diff_unpredict}], 1, [1 2], light_blue_color, blue_color)
axis([0.75 2.25 -2.5 7]); %axis square
plot(xlim, [0 0], 'k--')
xticks(1:2); xticklabels({'1', '6'}); xlabel('Problem')
ylabel('Change in rich-tone preference (d`)')
title('Amount of learning')
set(gcf, 'Position', [680   558   337   420])
%var_name = 'LearningAmount_compare'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')


    % stats: mixed model test for interaction
    % repeated measures anova
    datamtx = [[first_day_diff_predict last_day_diff_predict]; [first_day_diff_unpredict last_day_diff_unpredict]];
    between_subjs = [zeros(40,1); ones(40,1)];
    rmtbl = simple_mixed_anova(datamtx, between_subjs, {'problem'}, {'group'});
    ranova_amount_stats_comp = [rmtbl.DF(3) rmtbl.DF(5) rmtbl.F(5) rmtbl.pValue(5)];
    
    % ttest2s for group comparisons
    [~, pval_amount_comp_problem1, ~, stats_amount_comp_problem1] = ttest2(first_day_diff_predict, first_day_diff_unpredict);
    [~, pval_amount_comp_problem6, ~, stats_amount_comp_problem6] = ttest2(last_day_diff_predict, last_day_diff_unpredict);
    
    % compare to zero
    [~, pval_amount_predict_problem6_0, ~, stats_amount_predict_problem6_0] = ttest(last_day_diff_predict);
    [~, pval_amount_unpredict_problem6_0, ~, stats_amount_unpredict_problem6_0] = ttest(last_day_diff_unpredict);
    
    % ttests for compare problem 1 and problem 6 within group
    [~, pval_amount_predict_problems1and6, ~, stats_amount_predict_problems1and6] = ttest2(first_day_diff_predict, last_day_diff_predict)
    [~, pval_amount_unpredict_problems1and6, ~, stats_amount_unpredict_problems1and6] = ttest2(first_day_diff_unpredict, last_day_diff_unpredict)
    
    
   
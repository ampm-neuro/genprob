green_color = [21 154 72]./255; blue_color =  [43 73 159]./255;
light_green_color = [162 201 129]./255; light_blue_color = [102 164 213]./255;
figure_path = 'C:\Users\ampm1\Documents\manuscripts\gentrain\manuscript_figures\sfig02\';

%% Discrimination performance (see Figure1)
[firstLast_dprime_mtx_predict, ~, days_to_crit_predict] = ALL_stage_learn(['two_tone\' 'train_mevar_imaging_hpc'], 1:6, 2);
[firstLast_dprime_mtx_unpredict, ~, days_to_crit_unpredict] = ALL_stage_learn(['two_tone\' 'train_hivar_imaging_hpc'], 1:6, 2);
%for iplot = 1:14; close; end % close legacy plots

% predict
    figure; hold on
    for istage = 1:size(firstLast_dprime_mtx_predict,3)
        errorbar_plot([{firstLast_dprime_mtx_predict(:,1,istage)} {firstLast_dprime_mtx_predict(:,end,istage)}], 1, [-0.3 +0.3]+istage, light_green_color, green_color)
    end
    xlabel('Training problem (first and last day)')
    ylabel('Preference for rich tone (d`)')
    axis([0.25 6.75 -3 5]); axis square
    plot(xlim, [0 0], 'k--')
    title('Predictable learning curve')
    var_name = 'LearningCurve_predictable'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
    
        % stats: first days compared to day1
        % rm anova
        rmtbl_onewayanova_firstdaydiscrim_predict = simple_mixed_anova(squeeze(firstLast_dprime_mtx_predict(:,1,:)))
        
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
    axis([0.25 6.75 -3 5]); axis square
    plot(xlim, [0 0], 'k--')
	title('Predictable learning curve')
    var_name = 'LearningCurve_unpredictable'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
    
        % stats: first days compared to day1
        % repeated measures anova
        rmtbl_onewayanova_firstdaydiscrim_unpredict = simple_mixed_anova(squeeze(firstLast_dprime_mtx_unpredict(:,1,:)))
        
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

        
        
%% days to criterion line plots  (see Figure1)
        
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
    var_name = 'LearningRates_predictable'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
    
        % stats: first days compared to day1
        % repeated measures anova
        rmtbl_onewayanova_days2crit_predict = simple_mixed_anova(days_to_crit_predict)
        
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
    var_name = 'LearningRates_unpredictable'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
    
        % stats: first days compared to day1
        % repeated measures anova
        rmtbl_onewayanova_days2crit_unpredict = simple_mixed_anova(days_to_crit_unpredict)
        
        % pairwise ttests
        pvals_d2c_unpredict = nan(1,5); stats_cell_d2c_unpredict_d1 = cell(1,5);
        for itest = 1:5
            [~, pvals_d2c_unpredict(itest), ~, stats_cell_d2c_unpredict_d1{itest}] = ttest(days_to_crit_unpredict(:,1), days_to_crit_unpredict(:,itest+1));
        end
        
%% direct comparison plots  (see Figure1)

% first day
figure; hold on
errorbar_plot([{firstLast_dprime_mtx_predict(:,1,1)} {firstLast_dprime_mtx_predict(:,1,end)}], 1, [1 2], light_green_color, green_color)
errorbar_plot([{firstLast_dprime_mtx_unpredict(:,1,1)} {firstLast_dprime_mtx_unpredict(:,1,end)}], 1, [1 2], light_blue_color, blue_color)
axis([0.75 2.25 -2 3]); %axis square
yticks(-2:1:3)
plot(xlim, [0 0], 'k--')
xticks(1:2); xticklabels({'1', '6'}); xlabel('Problem')
ylabel('Preference for rich tone (d`)')
title('First day discrimination')
set(gcf, 'Position', [680   558   337   420])
var_name = 'LearningCurve_compare'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')

    %  model for interaction
    %
    datamtx = squeeze([firstLast_dprime_mtx_predict(:, 1, [1 6]); firstLast_dprime_mtx_unpredict(:, 1, [1 6])]);
    num_subjs_predict = 8; num_subjs_unpredict = 5; num_subjects = num_subjs_predict+num_subjs_unpredict; num_problems = 2;
    between_subjs = [zeros(num_subjs_predict,num_problems); ones(num_subjs_unpredict,num_problems)];
    subj_num_mtx = repmat((1:num_subjects)', 1, num_problems); % subject matrix
    stage_num_mtx = repmat(1:2, num_subjects, 1); % stage matrix
    tbl = table(datamtx(:), between_subjs(:), stage_num_mtx(:), subj_num_mtx(:),...
        'VariableNames',{'first_day_discrim', 'training_group', 'problem', 'subject'});
    tbl.subject = categorical(tbl.subject);
    model_str = 'first_day_discrim~training_group*problem+(1|subject)';
    lme_comp_firstdaydiscrim = fitlme(tbl, model_str)
    
    % ttest2s for group comparisons
    [~, pval_firstday_comp_problem1, ~, stats_firstday_comp_problem1] = ttest2(firstLast_dprime_mtx_predict(:,1,1), firstLast_dprime_mtx_unpredict(:,1,1))
    [~, pval_firstday_comp_problem6, ~, stats_firstday_comp_problem6] = ttest2(firstLast_dprime_mtx_predict(:,1,6), firstLast_dprime_mtx_unpredict(:,1,6))
    

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
var_name = 'LearningRate_compare'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')


    % model for interaction
    datamtx = [days_to_crit_predict(:, [1 6]); days_to_crit_unpredict(:, [1 6])];
    num_subjs_predict = 8; num_subjs_unpredict = 5; num_subjects = num_subjs_predict+num_subjs_unpredict; num_problems = 2;
    between_subjs = [zeros(num_subjs_predict,num_problems); ones(num_subjs_unpredict,num_problems)];
    subj_num_mtx = repmat((1:num_subjects)', 1, num_problems); % subject matrix
    stage_num_mtx = repmat(1:2, num_subjects, 1); % stage matrix
    tbl = table(datamtx(:), between_subjs(:), stage_num_mtx(:), subj_num_mtx(:),...
        'VariableNames',{'days2crit', 'training_group', 'problem', 'subject'});
    tbl.subject = categorical(tbl.subject);
    model_str = 'days2crit~training_group*problem+(1|subject)';
    lme_comp_days2crit = fitlme(tbl, model_str)
    

    % ttest2s for group comparisons
    [~, pval_d2c_comp_problem1, ~, stats_d2c_comp_problem1] = ttest2(days_to_crit_predict(:,1), days_to_crit_unpredict(:,1))
    [~, pval_d2c_comp_problem6, ~, stats_d2c_comp_problem6] = ttest2(days_to_crit_predict(:,6), days_to_crit_unpredict(:,6))
   

% amount of learning
figure; hold on
first_day_diff_predict = firstLast_dprime_mtx_predict(:,2,1) - firstLast_dprime_mtx_predict(:,1,1);
last_day_diff_predict = firstLast_dprime_mtx_predict(:,2,6) - firstLast_dprime_mtx_predict(:,1,6);
first_day_diff_unpredict = firstLast_dprime_mtx_unpredict(:,2,1) - firstLast_dprime_mtx_unpredict(:,1,1);
last_day_diff_unpredict = firstLast_dprime_mtx_unpredict(:,2,6) - firstLast_dprime_mtx_unpredict(:,1,6);
errorbar_plot([{first_day_diff_predict} {last_day_diff_predict}], 1, [1 2], light_green_color, green_color)
errorbar_plot([{first_day_diff_unpredict} {last_day_diff_unpredict}], 1, [1 2], light_blue_color, blue_color)
axis([0.75 2.25 -2 5]); %axis square
plot(xlim, [0 0], 'k--')
xticks(1:2); xticklabels({'1', '6'}); xlabel('Problem')
ylabel('Change in rich-tone preference (d`)')
title('Amount of learning')
set(gcf, 'Position', [680   558   337   420])
var_name = 'LearningAmount_compare'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')


    % stats: test for interaction
    % repeated measures anova
    datamtx = [[first_day_diff_predict last_day_diff_predict]; [first_day_diff_unpredict last_day_diff_unpredict]];
    num_subjs_predict = 8; num_subjs_unpredict = 5; num_subjects = num_subjs_predict+num_subjs_unpredict; num_problems = 2;
    between_subjs = [zeros(num_subjs_predict,num_problems); ones(num_subjs_unpredict,num_problems)];
    subj_num_mtx = repmat((1:num_subjects)', 1, num_problems); % subject matrix
    stage_num_mtx = repmat(1:2, num_subjects, 1); % stage matrix
    tbl = table(datamtx(:), between_subjs(:), stage_num_mtx(:), subj_num_mtx(:),...
        'VariableNames',{'amount_learn', 'training_group', 'problem', 'subject'});
    tbl.subject = categorical(tbl.subject);
    model_str = 'amount_learn~training_group*problem+(1|subject)';
    lme_comp_amountlearn = fitlme(tbl, model_str);
    
    % ttest2s for group comparisons
    [~, pval_amount_comp_problem1, ~, stats_amount_comp_problem1] = ttest2(first_day_diff_predict, first_day_diff_unpredict);
    [~, pval_amount_comp_problem6, ~, stats_amount_comp_problem6] = ttest2(last_day_diff_predict, last_day_diff_unpredict);
    
    % compare unpredict to zero
    [~, pval_amount_predict_problem6_0, ~, stats_amount_predict_problem6_0] = ttest(last_day_diff_predict);
    
    % ttests for compare problem 1 and problem 6 within group
    [~, pval_amount_predict_problems1and6, ~, stats_amount_predict_problems1and6] = ttest2(first_day_diff_predict, last_day_diff_predict);
    [~, pval_amount_unpredict_problems1and6, ~, stats_amount_unpredict_problems1and6] = ttest2(first_day_diff_unpredict, last_day_diff_unpredict);
    
    
    
  
    
%% all probe curves (see Figure2)
figure;

% predict
probe_fps_predict = [get_file_paths_targeted('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_imaging_hpc', 'probe')];            
for iprobe = 1:7
    
    if iprobe<7
        probe_paths = probe_fps_predict(contains(probe_fps_predict, ['preprobe_0' num2str(iprobe)]));
    else
        probe_paths = probe_fps_predict(contains(probe_fps_predict, 'postprobe_01'));
    end
    subplot(2,7,iprobe)
    plot_allprobes(probe_paths); 
    rich_bounds_prob('mevar', iprobe-1);
    ylim([0 60]); xlim([4500 38500])
    xticks([]); xticklabels([])
    xlabel([])
    if iprobe==1
        yticks([0 60]);
    else
        yticks([]); yticklabels([])
        ylabel([])
    end
    axis square
    title(['Probe ' num2str(iprobe-1)])
end

% unpredict
probe_fps_predict = [get_file_paths_targeted('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_hivar_imaging_hpc', 'probe')];            
for iprobe = 1:7
    
    if iprobe<7
        probe_paths = probe_fps_predict(contains(probe_fps_predict, ['preprobe_0' num2str(iprobe)]));
    else
        probe_paths = probe_fps_predict(contains(probe_fps_predict, 'postprobe_01'));
    end
    subplot(2,7,iprobe+7)
    plot_allprobes(probe_paths); 
    rich_bounds_prob('hivar', iprobe-1);
    ylim([0 60]); xlim([4500 38500])
    yticks([0 60]); xticks([5000 35000]); xticklabels({'5', '35'})
    xlabel('Tone freq. (kHz)')
    axis square
    if iprobe==1
        yticks([0 60]);
    else
        yticks([]); yticklabels([])
        ylabel([])
    end
    title([])
end

% save
set(gcf, 'Position', [10 501 1905 491])
var_name = 'probe_curves'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')



%% correlation matrices and associated line plots  (see Figure2)
[probe_waits_predict, ~, cm_17_predict, offdiag_cell_predict, comp2one_cell_predict, comp2seven_cell_predict] = probe_corr_mtx('train_mevar_imaging_hpc');
%var_name = 'ProbeBeh_corm_predict'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
[probe_waits_unpredict, ~, cm_17_unpredict, offdiag_cell_unpredict, comp2one_cell_unpredict, comp2seven_cell_unpredict] = probe_corr_mtx('train_hivar_imaging_hpc');
%var_name = 'ProbeBeh_corm_unpredict'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')

% compare similarity of responses on adjacent probes
% predict
figure; hold on; 
errorbar_plot(offdiag_cell_unpredict, 1, [], light_blue_color, blue_color);
plot(xlim, [1 1].*0, 'k--')
xticklabels({'0&1', '1&2', '2&3', '3&4', '4&5', '5&6'}); 
xlabel('Probe comparison')
ylabel('Similarity (r)')
ylim([-1 1])
title('Compare adjacent probes')
axis square
var_name = 'ProbeBeh_adjacent_predict'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')

        % stats: first days compared to day1
        % repeated measures anova
        rmtbl = simple_mixed_anova(cell2mat(offdiag_cell_predict));
        ranova_adjacent_stats_predict = [rmtbl.DF(2) rmtbl.DF(3) rmtbl.F(3) rmtbl.pValue(3)];
        
        % pairwise ttests
        pvals_adjacent_predict_d1 = nan(1,5); stats_cell_adjacent_predict_d1 = cell(1,5);
        for itest = 1:5
            [~, pvals_adjacent_predict_d1(itest), ~, stats_cell_adjacent_predict_d1{itest}] = ttest(offdiag_cell_predict{1}, offdiag_cell_predict{itest+1});
        end
        
% unpredict
figure; hold on; 
errorbar_plot(offdiag_cell_predict, 1, [], light_green_color, green_color);
plot(xlim, [1 1].*0, 'k--')
xticklabels({'0&1', '1&2', '2&3', '3&4', '4&5', '5&6'}); 
xlabel('Probe comparison')
ylabel('Similarity (r)')
ylim([-1 1])
title('Compare adjacent probes')
axis square
var_name = 'ProbeBeh_adjacent_unpredict'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')

        % stats: first days compared to day1
        % repeated measures anova
        rmtbl = simple_mixed_anova(cell2mat(offdiag_cell_unpredict));
        ranova_adjacent_stats_unpredict = [rmtbl.DF(2) rmtbl.DF(3) rmtbl.F(3) rmtbl.pValue(3)];
        
        % pairwise ttests
        pvals_adjacent_unpredict_d1 = nan(1,5); stats_cell_adjacent_unpredict_d1 = cell(1,5);
        for itest = 1:5
            [~, pvals_adjacent_unpredict_d1(itest), ~, stats_cell_adjacent_unpredict_d1{itest}] = ttest(offdiag_cell_unpredict{1}, offdiag_cell_unpredict{itest+1});
        end
        
% interaction compare first and last
figure; hold on; 
errorbar_plot(offdiag_cell_unpredict([1 end]), 1, [], light_blue_color, blue_color);
errorbar_plot(offdiag_cell_predict([1 end]), 1, [], light_green_color, green_color);
plot(xlim, [1 1].*0, 'k--')
xticks(1:2); xticklabels({'0&1', '5&6'});
xlabel('Probe comparison')
ylabel('Similarity (r)')
ylim([-1 1]); xlim([0.75 2.25])
set(gcf, 'Position', [680   558   337   420])
var_name = 'ProbeBeh_adjacent_comp'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')

    %  stats: mixed model test for interaction
    % repeated measures anova
    datamtx = [cell2mat(offdiag_cell_predict([1 end])); cell2mat(offdiag_cell_unpredict([1 end]))];
    between_subjs = [zeros(40,1); ones(40,1)];
    rmtbl = simple_mixed_anova(datamtx, between_subjs, {'problem'}, {'group'});
    ranova_amount_stats_comp = [rmtbl.DF(3) rmtbl.DF(5) rmtbl.F(5) rmtbl.pValue(5)];
    
    % ttest2s for group comparisons
    [~, pval_amount_comp_problem1, ~, stats_amount_comp_problem1] = ttest2(cell2mat(offdiag_cell_predict(1)), cell2mat(offdiag_cell_unpredict(1)));
    [~, pval_amount_comp_problem6, ~, stats_amount_comp_problem6] = ttest2(cell2mat(offdiag_cell_predict(end)), cell2mat(offdiag_cell_unpredict(end)));
    
    % compare unpredict to zero
    [~, pval_amount_predict_problem1_0, ~, stats_amount_predict_problem1_0] = ttest(cell2mat(offdiag_cell_predict(1)));
    [~, pval_amount_predict_problem6_0, ~, stats_amount_predict_problem6_0] = ttest(cell2mat(offdiag_cell_predict(end)));
    [~, pval_amount_unpredict_problem1_0, ~, stats_amount_unpredict_problem1_0] = ttest(cell2mat(offdiag_cell_unpredict(1)));
    [~, pval_amount_unpredict_problem6_0, ~, stats_amount_unpredict_problem6_0] = ttest(cell2mat(offdiag_cell_unpredict(end)));
    
    

%% preference for rich tones data (see fig3)
%
[~, ~, diff_waits_problem_preprobe_predict, diff_waits_problem_postprobe_predict, correlation_plot_cell_predict] = ALL_probe_vs_behavior_fit('train_mevar_imaging_hpc');
%for ifig = 1:34; close; end
[~, ~, diff_waits_problem_preprobe_unpredict, diff_waits_problem_postprobe_unpredict, correlation_plot_cell_unpredict] = ALL_probe_vs_behavior_fit('train_hivar_imaging_hpc');
%for ifig = 1:34; close; end

    previous_predict = cell(1,size(diff_waits_problem_postprobe_predict,2));
    previous_unpredict = cell(1,size(diff_waits_problem_postprobe_unpredict,2));
    next_predict = cell(1,size(diff_waits_problem_preprobe_predict,2));
    next_unpredict = cell(1,size(diff_waits_problem_preprobe_unpredict,2));
    for ic = 1:size(diff_waits_problem_postprobe_predict,2)
        previous_predict{ic} = diff_waits_problem_postprobe_predict(:,ic);
        previous_unpredict{ic} = diff_waits_problem_postprobe_unpredict(:,ic);
        next_predict{ic} = diff_waits_problem_preprobe_predict(:,ic);
        next_unpredict{ic} = diff_waits_problem_preprobe_unpredict(:,ic);
    end
save('sfigure02_2', 'previous_predict', 'previous_unpredict', 'next_predict', 'next_unpredict', 'diff_waits_problem_preprobe_predict', 'diff_waits_problem_postprobe_predict', 'diff_waits_problem_preprobe_unpredict', 'diff_waits_problem_postprobe_unpredict', 'correlation_plot_cell_predict', 'correlation_plot_cell_unpredict')
%}
load('sfigure02_2', 'previous_predict', 'previous_unpredict', 'next_predict', 'next_unpredict', 'diff_waits_problem_preprobe_predict', 'diff_waits_problem_postprobe_predict', 'diff_waits_problem_preprobe_unpredict', 'diff_waits_problem_postprobe_unpredict', 'correlation_plot_cell_predict', 'correlation_plot_cell_unpredict')

%% preference for rich tones plot (see fig3)
% previous problem
    
    % predict
    figure; hold on
    errorbar_plot(previous_predict, 1, 2:7, light_green_color, green_color)
    xlim([1.5 7.5]); ylim([-20 25]); axis square
    xlabel('Probe'); xticks(2:7); xticklabels(1:6)
    ylabel('Rich - poor tone wait (s)')
    plot(xlim, [1 1].*0, 'k--')
    title('Previous problem')
    var_name = 'ProbeBeh_ToneDiscr_PastProblem_predict'; 
    print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
    
        % ttest compared to zero
        pvals_previous_predict_0 = nan(size(previous_predict,2),1);
        stats_previous_predict_0 = cell(size(pvals_previous_predict_0));
        for itest = 1:size(previous_predict,2)
            [~, pvals_previous_predict_0(itest), ~, stats_previous_predict_0{itest}] = ttest(cell2mat(previous_predict(itest)));
        end
    
    % unpredict
    figure; hold on
    errorbar_plot(previous_unpredict, 1, 2:7, light_blue_color, blue_color)
    xlim([1.5 7.5]); ylim([-20 25]); axis square
    xlabel('Probe'); xticks(2:7); xticklabels(1:6)
    ylabel('Rich - poor tone wait (s)')
    plot(xlim, [1 1].*0, 'k--')
    title('Previous problem')
    var_name = 'ProbeBeh_ToneDiscr_PastProblem_unpredict'; 
    print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
    
        % ttest compared to zero
        pvals_previous_unpredict_0 = nan(size(previous_unpredict,2),1);
        stats_previous_unpredict_0 = cell(size(pvals_previous_unpredict_0));
        for itest = 1:size(previous_predict,2)
            [~, pvals_previous_unpredict_0(itest), ~, stats_previous_unpredict_0{itest}] = ttest(cell2mat(previous_unpredict(itest)));
        end
        
    % compare
    figure; hold on
    errorbar_plot(previous_unpredict([1 end]), 1, 1:2, light_blue_color, blue_color)
    errorbar_plot(previous_predict([1 end]), 1, 1:2, light_green_color, green_color)
    xlim([0.75 2.25]); ylim([-20 25]); set(gcf, 'Position', [680   558   337   420])
    xticks(1:2); xticklabels({'1', '6'})
    plot(xlim, [1 1].*0, 'k--')
    xlabel('Probe')
    ylabel('Rich - poor tone wait (s)')
    title('Previous problem')
    var_name = 'ProbeBeh_ToneDiscr_PastProblem_comp'; 
    print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
    
        %  stats: mixed model test for interaction
        % repeated measures anova
        datamtx = [cell2mat(previous_unpredict([1 end])); cell2mat(previous_predict([1 end]))];
        between_subjs = [zeros(length(previous_unpredict{1}),1); ones(length(previous_predict{1}),1)];
        rmtbl = simple_mixed_anova(datamtx, between_subjs, {'problem'}, {'group'});
        ranova_previous_stats_comp = [rmtbl.DF(3) rmtbl.DF(5) rmtbl.F(5) rmtbl.pValue(5)];

        % ttest2s for group comparisons
        [~, pval_previous_comp_problem1, ~, stats_amount_comp_problem1] = ttest2(cell2mat(previous_predict(1)), cell2mat(previous_unpredict(1)));
        [~, pval_previous_comp_problem6, ~, stats_amount_comp_problem6] = ttest2(cell2mat(previous_predict(end)), cell2mat(previous_unpredict(end)));
    
    
    
%next problem

    % predict
    figure; hold on
    errorbar_plot(next_predict, 1, 1:6, light_green_color, green_color)
    xlim([0.5 6.5]); ylim([-20 25]); axis square
    xlabel('Probe'); xticks(1:6); xticklabels(0:5)
    ylabel('Rich - poor tone wait (s)')
    plot(xlim, [1 1].*0, 'k--')
    title('Next problem')
    var_name = 'ProbeBeh_ToneDiscr_NextProblem_predict'; 
    print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
        
        % ttest compared to zero
        pvals_next_predict_0 = nan(size(next_predict,2),1);
        stats_next_predict_0 = cell(size(pvals_next_predict_0));
        for itest = 1:size(next_predict,2)
            [~, pvals_next_predict_0(itest), ~, stats_next_predict_0{itest}] = ttest(cell2mat(next_predict(itest)));
        end
        
    % unpredict
    figure; hold on
    errorbar_plot(next_unpredict, 1, 1:6, light_blue_color, blue_color)
    xlim([0.5 6.5]); ylim([-20 25]); axis square
    xlabel('Probe'); xticks(1:6); xticklabels(0:5)
    ylabel('Rich - poor tone wait (s)')
    plot(xlim, [1 1].*0, 'k--')
    title('Next problem')
    var_name = 'ProbeBeh_ToneDiscr_NextProblem_unpredict'; 
    print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
    
        % ttest compared to zero
        pvals_next_unpredict_0 = nan(size(next_unpredict,2),1);
        stats_next_unpredict_0 = cell(size(pvals_next_unpredict_0));
        for itest = 1:size(next_unpredict,2)
            [~, pvals_next_unpredict_0(itest), ~, stats_next_unpredict_0{itest}] = ttest(cell2mat(next_unpredict(itest)));
        end
    
    % compare
    figure; hold on
    errorbar_plot(next_unpredict([1 end]), 1, 1:2, light_blue_color, blue_color)
    errorbar_plot(next_predict([1 end]), 1, 1:2, light_green_color, green_color)
    xlim([0.75 2.25]); ylim([-20 25]); set(gcf, 'Position', [680   558   337   420])
    xticks(1:2); xticklabels({'0', '5'})
    plot(xlim, [1 1].*0, 'k--')
    xlabel('Probe')
    ylabel('Rich - poor tone wait (s)')
    title('Next problem')
    var_name = 'ProbeBeh_ToneDiscr_NextProblem_comp'; 
    print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
    
        %  stats: mixed model test for interaction
        % repeated measures anova
        datamtx = [cell2mat(next_unpredict([1 end])); cell2mat(next_predict([1 end]))];
        between_subjs = [zeros(length(next_unpredict{1}),1); ones(length(next_predict{1}),1)];
        rmtbl = simple_mixed_anova(datamtx, between_subjs, {'problem'}, {'group'});
        ranova_next_stats_comp = [rmtbl.DF(3) rmtbl.DF(5) rmtbl.F(5) rmtbl.pValue(5)];

        % ttest2s for group comparisons
        [~, pval_next_comp_problem1, ~, stats_next_comp_problem1] = ttest2(cell2mat(next_predict(1)), cell2mat(next_unpredict(1)));
        [~, pval_next_comp_problem6, ~, stats_next_comp_problem6] = ttest2(cell2mat(next_predict(end)), cell2mat(next_unpredict(end)));
    
    
    
    
%% correlation plots

% combine across groups
accuracy_of_prediction = [diff_waits_problem_preprobe_predict; diff_waits_problem_preprobe_unpredict];
first_day_discrimination = [correlation_plot_cell_predict{1}; correlation_plot_cell_unpredict{1}];
days_to_criterion = [correlation_plot_cell_unpredict{2}; correlation_plot_cell_predict{2}];
amount_learning = [correlation_plot_cell_predict{3}; correlation_plot_cell_unpredict{3}];
probe_over_probe_corr = [correlation_plot_cell_predict{4}; correlation_plot_cell_unpredict{4}];

    % first day discrimination
    plot2var(diff_waits_problem_preprobe_predict, correlation_plot_cell_predict{1}, ...
            diff_waits_problem_preprobe_unpredict, correlation_plot_cell_unpredict{1}, 'Prediction accuracy and initial discrimination')
        axis([-25 25 -3 5]); ylabel('Discrimination (d`)')
        hold on; plot(xlim, [0 0], 'k--'); plot([0 0], ylim, 'k--')
        var_name = 'PredictAccuracy_Firstday_corr'; 
    	print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
    % last day discrimination
    plot2var(diff_waits_problem_preprobe_predict, correlation_plot_cell_predict{5}, ...
            diff_waits_problem_preprobe_unpredict, correlation_plot_cell_unpredict{5}, 'Prediction accuracy and final discrimination')
        axis([-25 25 -3 5]); ylabel('Discrimination (d`)')
        hold on; plot(xlim, [0 0], 'k--'); plot([0 0], ylim, 'k--')
        var_name = 'PredictAccuracy_Lastday_corr'; 
    	print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
    % days to criterion 
    plot2var(diff_waits_problem_preprobe_predict, correlation_plot_cell_predict{2}, ...
            diff_waits_problem_preprobe_unpredict, correlation_plot_cell_unpredict{2}, 'Prediction accuracy and learning rate')
        axis([-25 25 .5 10.5]); ylabel('Days to criterion')
        hold on; plot(xlim, [0 0], 'k--'); plot([0 0], ylim, 'k--')
        var_name = 'PredictAccuracy_Days2crit_corr'; 
    	print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
    % amount learning 
    plot2var(diff_waits_problem_preprobe_predict, correlation_plot_cell_predict{3}, ...
            diff_waits_problem_preprobe_unpredict, correlation_plot_cell_unpredict{3}, 'Prediction accuracy and amount of learning')
        axis([-25 25 -2 5]); ylabel('Last day - first day (d`)')
        hold on; plot(xlim, [0 0], 'k--'); plot([0 0], ylim, 'k--')
        var_name = 'PredictAccuracy_AmountLearn_corr'; 
    	print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
    % memory stability
    plot2var(diff_waits_problem_preprobe_predict, correlation_plot_cell_predict{4}, ...
            diff_waits_problem_preprobe_unpredict, correlation_plot_cell_unpredict{4}, 'Prediction accuracy and memory stability')
        axis([-25 25 -.6 .6]); ylabel('Correlation between pre and post probe (r)')
        hold on; plot(xlim, [0 0], 'k--'); plot([0 0], ylim, 'k--')
        var_name = 'PredictAccuracy_MemoryStability_corr'; 
    	print([figure_path var_name], '-dpdf', '-painters', '-bestfit')

%mixed model prep
num_subjects = 13; num_problems = 6;
subj_num_mtx = repmat((1:num_subjects)', 1, num_problems); % subject matrix
stage_num_mtx = repmat(1:6, num_subjects, 1); % stage matrix
model_str = 'var1~var2+(1|subject)';

% mixed models
lme_first_day_discrim = lme_function(model_str, first_day_discrimination, accuracy_of_prediction, subj_num_mtx, stage_num_mtx)
lme_days2crit = lme_function(model_str, days_to_criterion, accuracy_of_prediction, subj_num_mtx, stage_num_mtx)
lme_amount = lme_function(model_str, amount_learning, accuracy_of_prediction, subj_num_mtx, stage_num_mtx)
lme_stability = lme_function(model_str, probe_over_probe_corr, accuracy_of_prediction, subj_num_mtx, stage_num_mtx)

function lme = lme_function(model_str, var1, var2, subj_num_mtx, stage_num_mtx)
    % function testing whether var2 (input 3) predicts var1 (input 2)

    % update model string with input names

        %var1
        var1_str_starts = strfind(model_str, 'var1');
        for v1ss = fliplr(var1_str_starts)
            model_str = [model_str(1:v1ss-1) inputname(2) model_str((v1ss+length('var1')):end)];
        end

        %var2
        var2_str_starts = strfind(model_str, 'var2');
        for v2ss = fliplr(var2_str_starts)
            model_str = [model_str(1:v2ss-1) inputname(3) model_str((v2ss+length('var2')):end)];
        end
        
    % mixed model
    tbl = table(var1(:), var2(:), subj_num_mtx(:), stage_num_mtx(:),...
        'VariableNames',{inputname(2), inputname(3), 'subject', 'problem'});
    tbl.subject = categorical(tbl.subject);
    %tbl.problem = categorical(tbl.problem);
    lme = fitlme(tbl, model_str);
                
end

% Figure functions
function plot2var(var1_predict, var2_predict, var1_unpredict, var2_unpredict, title_str)
% vars need to be matrices with (subjects, samples)

    % colormap
    green_colors = [162 221 154; 118 205 117; 065 181 093; 037 149 069; 007 120 056; 013 078 032]./255;
    blue_colors = [158 202 225; 106 173 214; 067 146 198; 035 114 181; 009 083 157; 033 052 104]./255;
    
    % open figure
    figure; hold on
    
    % dot plot 
        %unpredict
        for isubj = 1:size(var1_unpredict,1)
            for iprobe = 1:size(var1_unpredict,2)
                plot(var1_unpredict(isubj, iprobe), var2_unpredict(isubj, iprobe), 'o', 'color', blue_colors(iprobe,:))
            end
        end
        
        % predict
        for isubj = 1:size(var1_predict,1)
            for iprobe = 1:size(var1_predict,2)
                plot(var1_predict(isubj, iprobe), var2_predict(isubj, iprobe), 'o', 'color', green_colors(iprobe,:))
            end
        end
        
    % means 
        % unpredict
        for iprobe = 1:size(var1_unpredict,2)
            plot(nanmean(var1_unpredict(:, iprobe)), nanmean(var2_unpredict(:, iprobe)), '.', 'color', blue_colors(iprobe,:), 'markersize', 60)
        end
        % predict
        for iprobe = 1:size(var1_predict,2)
            plot(nanmean(var1_predict(:, iprobe)), nanmean(var2_predict(:, iprobe)), '.', 'color', green_colors(iprobe,:), 'markersize', 60)
        end
    
    [r, p] = fit_line([var1_predict(:); var1_unpredict(:)], [var2_predict(:); var2_unpredict(:)], 0);
    title([title_str '; r=' num2str(r) ', p=' num2str(p)])
    set(gca,'TickLength',[0, 0]); box off;
    xlabel('Future rich tone wait - future poor tone wait (s)')
    axis square
    

end


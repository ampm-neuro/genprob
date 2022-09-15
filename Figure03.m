green_color = [21 154 72]./255; blue_color =  [43 73 159]./255;
light_green_color = [162 201 129]./255; light_blue_color = [102 164 213]./255;
figure_path = 'C:\Users\ampm1\Documents\manuscripts\gentrain\manuscript_figures\fig03\';


%% memory for recent discrimination example diagram

% search for examples
%{
data_path_mevar = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar';
subject_folders = get_folder_paths_all(data_path_mevar);
% iterate through subjects
for isubj = 1:size(subject_folders,1)
    figure
    
    % load last problem 1 session
    problem1 = get_file_paths_targeted(subject_folders{isubj}, 'gen07');
    load(problem1{end}, 'trl_mtx')
    
    % plot 
    try
        subplot(1,2,1); hold on
        wait_times_plot(trl_mtx, 3);
        ylim([0 60])
    catch; disp('catch'); end
    
    % load last probe 1 session
    probe1 = get_file_paths_targeted(subject_folders{isubj}, 'preprobe_02');
    load(probe1{1}, 'trl_mtx')
    
    % plot if 82 trials
    if size(trl_mtx,1)>=82
        try
            subplot(1,2,2); hold on
            wait_times_plot(trl_mtx); wait_times_plot_tonecolor(trl_mtx);
            ylim([0 60])
            sgtitle(num2str(isubj));
            drawnow
        catch; disp('catch'); close; end
    end
end
%}

%
figure;
load('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar\687930m1\gen07_mevar01-10.mat', 'trl_mtx')
subplot(1,2,1); wait_times_plot(trl_mtx, 3); 
ylim([0 40]); rich_bounds_prob('mevar', 1); xlim([4500 38500]); axis square
load('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar\687930m1\gen14_preprobe_02_03d.mat', 'trl_mtx')
subplot(1,2,2); hold on; wait_times_plot(trl_mtx); wait_times_plot_tonecolor(trl_mtx); 
ylim([0 40]); rich_bounds_prob('mevar', 1); xlim([4500 38500]); axis square
sgtitle('687930m1');
var_name = 'ProbeBeh_eg_recentdiscrim'; 
%print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
%}



%% memory for future discrimination example diagram

% search for examples
%{
data_path_mevar = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar';
subject_folders = get_folder_paths_all(data_path_mevar);
% iterate through subjects
for isubj = 1:size(subject_folders,1)
    figure
    
    % load last problem 1 session
    problem6 = get_file_paths_targeted(subject_folders{isubj}, 'gen12');
    load(problem6{1}, 'trl_mtx')
    
    % plot 
    try
        subplot(1,2,2); hold on
        wait_times_plot(trl_mtx, 3);
        ylim([0 60])
    catch; disp('catch'); end
    
    % load last probe 5 session
    probe1 = get_file_paths_targeted(subject_folders{isubj}, 'preprobe_06');
    load(probe1{1}, 'trl_mtx')
    
    % plot if 82 trials
    if size(trl_mtx,1)>=82
        try
            subplot(1,2,1); hold on
            wait_times_plot(trl_mtx); wait_times_plot_tonecolor(trl_mtx);
            ylim([0 60])
            sgtitle(num2str(isubj));
            drawnow
        catch; disp('catch'); close; end
    end
end
%}

%
figure;
load('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar\682611m3\gen14_preprobe_06_01d.mat', 'trl_mtx')
subplot(1,2,1); hold on; wait_times_plot(trl_mtx); wait_times_plot_tonecolor(trl_mtx); 
ylim([0 40]); rich_bounds_prob('mevar', 6); xlim([4500 38500]); axis square
load('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar\682611m3\gen12_mevar06-01.mat', 'trl_mtx')
subplot(1,2,2); wait_times_plot(trl_mtx, 3); 
ylim([0 40]); rich_bounds_prob('mevar', 6); xlim([4500 38500]); axis square
sgtitle('682611m3');
var_name = 'ProbeBeh_eg_futurediscrim'; 
%print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
%}



%% preference for rich tones data
%{
[~, ~, diff_waits_problem_preprobe_predict, diff_waits_problem_postprobe_predict, correlation_plot_cell_predict] = ALL_probe_vs_behavior_fit('train_mevar');
%for ifig = 1:34; close; end
[~, ~, diff_waits_problem_preprobe_unpredict, diff_waits_problem_postprobe_unpredict, correlation_plot_cell_unpredict] = ALL_probe_vs_behavior_fit('train_hivar');
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
    save('Figure04', 'previous_predict', 'previous_unpredict', 'next_predict', 'next_unpredict', 'diff_waits_problem_preprobe_predict', 'diff_waits_problem_postprobe_predict', 'diff_waits_problem_preprobe_unpredict', 'diff_waits_problem_postprobe_unpredict', 'correlation_plot_cell_predict', 'correlation_plot_cell_unpredict')
%}
load('Figure04', 'previous_predict', 'previous_unpredict', 'next_predict', 'next_unpredict', 'diff_waits_problem_preprobe_predict', 'diff_waits_problem_postprobe_predict', 'diff_waits_problem_preprobe_unpredict', 'diff_waits_problem_postprobe_unpredict', 'correlation_plot_cell_predict', 'correlation_plot_cell_unpredict')
    

%% preference for rich tones plot
% previous problem
    
    % predict
    figure; hold on
    errorbar_plot(previous_predict, 1, 2:7, light_green_color, green_color)
    xlim([1.5 7.5]); ylim([-30 30]); axis square
    xlabel('Probe'); xticks(2:7); xticklabels(1:6)
    ylabel('Rich - poor tone wait (s)')
    plot(xlim, [1 1].*0, 'k--')
    title('Previous problem')
    var_name = 'ProbeBeh_ToneDiscr_PastProblem_predict'; 
    %print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
    
        % ttest compared to zero
        pvals_previous_predict_0 = nan(size(previous_predict,2),1);
        stats_previous_predict_0 = cell(size(pvals_previous_predict_0));
        for itest = 1:size(previous_predict,2)
            [~, pvals_previous_predict_0(itest), ~, stats_previous_predict_0{itest}] = ttest(cell2mat(previous_predict(itest)));
        end
    
    % unpredict
    figure; hold on
    errorbar_plot(previous_unpredict, 1, 2:7, light_blue_color, blue_color)
    xlim([1.5 7.5]); ylim([-30 30]); axis square
    xlabel('Probe'); xticks(2:7); xticklabels(1:6)
    ylabel('Rich - poor tone wait (s)')
    plot(xlim, [1 1].*0, 'k--')
    title('Previous problem')
    var_name = 'ProbeBeh_ToneDiscr_PastProblem_unpredict'; 
    %print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
    
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
    xlim([0.75 2.25]); ylim([-30 30]); set(gcf, 'Position', [680   558   337   420])
    xticks(1:2); xticklabels({'1', '6'})
    plot(xlim, [1 1].*0, 'k--')
    xlabel('Probe')
    ylabel('Rich - poor tone wait (s)')
    title('Previous problem')
    var_name = 'ProbeBeh_ToneDiscr_PastProblem_comp'; 
    %print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
    
        %  stats: mixed model test for interaction
        % repeated measures anova
        datamtx = [cell2mat(previous_unpredict([1 end])); cell2mat(previous_predict([1 end]))];
        between_subjs = [zeros(40,1); ones(40,1)];
        rmtbl = simple_mixed_anova(datamtx, between_subjs, {'problem'}, {'group'});
        ranova_previous_stats_comp = [rmtbl.DF(3) rmtbl.DF(5) rmtbl.F(5) rmtbl.pValue(5)];

        % ttest2s for group comparisons
        [~, pval_previous_comp_problem1, ~, stats_amount_comp_problem1] = ttest2(cell2mat(previous_predict(1)), cell2mat(previous_unpredict(1)));
        [~, pval_previous_comp_problem6, ~, stats_amount_comp_problem6] = ttest2(cell2mat(previous_predict(end)), cell2mat(previous_unpredict(end)));
    
    
    
%next problem

    % predict
    figure; hold on
    errorbar_plot(next_predict, 1, 1:6, light_green_color, green_color)
    xlim([0.5 6.5]); ylim([-30 30]); axis square
    xlabel('Probe'); xticks(1:6); xticklabels(0:5)
    ylabel('Rich - poor tone wait (s)')
    plot(xlim, [1 1].*0, 'k--')
    title('Next problem')
    var_name = 'ProbeBeh_ToneDiscr_NextProblem_predict'; 
    %print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
        
        % ttest compared to zero
        pvals_next_predict_0 = nan(size(next_predict,2),1);
        stats_next_predict_0 = cell(size(pvals_next_predict_0));
        for itest = 1:size(next_predict,2)
            [~, pvals_next_predict_0(itest), ~, stats_next_predict_0{itest}] = ttest(cell2mat(next_predict(itest)));
        end
        
    % unpredict
    figure; hold on
    errorbar_plot(next_unpredict, 1, 1:6, light_blue_color, blue_color)
    xlim([0.5 6.5]); ylim([-30 30]); axis square
    xlabel('Probe'); xticks(1:6); xticklabels(0:5)
    ylabel('Rich - poor tone wait (s)')
    plot(xlim, [1 1].*0, 'k--')
    title('Next problem')
    var_name = 'ProbeBeh_ToneDiscr_NextProblem_unpredict'; 
    %print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
    
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
    xlim([0.75 2.25]); ylim([-30 30]); set(gcf, 'Position', [680   558   337   420])
    xticks(1:2); xticklabels({'0', '5'})
    plot(xlim, [1 1].*0, 'k--')
    xlabel('Probe')
    ylabel('Rich - poor tone wait (s)')
    title('Next problem')
    var_name = 'ProbeBeh_ToneDiscr_NextProblem_comp'; 
    %print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
    
        %  stats: mixed model test for interaction
        % repeated measures anova
        datamtx = [cell2mat(next_unpredict([1 end])); cell2mat(next_predict([1 end]))];
        between_subjs = [zeros(40,1); ones(40,1)];
        rmtbl = simple_mixed_anova(datamtx, between_subjs, {'problem'}, {'group'});
        ranova_previous_stats_comp = [rmtbl.DF(3) rmtbl.DF(5) rmtbl.F(5) rmtbl.pValue(5)];

        % ttest2s for group comparisons
        [~, pval_next_comp_problem1, ~, stats_next_comp_problem1] = ttest2(cell2mat(next_predict(1)), cell2mat(next_unpredict(1)));
        [~, pval_next_comp_problem6, ~, stats_next_comp_problem6] = ttest2(cell2mat(next_predict(end)), cell2mat(next_unpredict(end)));
    
    
    
    
%% correlation plots

% combine across groups
accuracy_of_prediction = [diff_waits_problem_preprobe_predict; diff_waits_problem_preprobe_unpredict];
first_day_discrimination = [correlation_plot_cell_predict{1}; correlation_plot_cell_unpredict{1}];
days_to_criterion = [correlation_plot_cell_predict{2}; correlation_plot_cell_unpredict{2}];
amount_learning = [correlation_plot_cell_predict{3}; correlation_plot_cell_unpredict{3}];
probe_over_probe_corr = [correlation_plot_cell_predict{4}; correlation_plot_cell_unpredict{4}];
last_day_discrimination = [correlation_plot_cell_predict{5}; correlation_plot_cell_unpredict{5}];

    % first day discrimination
    plot2var(diff_waits_problem_preprobe_predict, correlation_plot_cell_predict{1}, ...
            diff_waits_problem_preprobe_unpredict, correlation_plot_cell_unpredict{1}, 'Prediction accuracy and initial discrimination')
        axis([-25 25 -3 5.5]); ylabel('First day discrimination (d`)')
        var_name = 'PredictAccuracy_Firstday_corr'; 
    	%print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
    plot2var(diff_waits_problem_preprobe_predict, correlation_plot_cell_predict{2}, ...
            diff_waits_problem_preprobe_unpredict, correlation_plot_cell_unpredict{2}, 'Prediction accuracy and learning rate')
        axis([-25 25 .5 10.5]); ylabel('Days to criterion')
        var_name = 'PredictAccuracy_Days2crit_corr'; 
    	%print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
    plot2var(diff_waits_problem_preprobe_predict, correlation_plot_cell_predict{3}, ...
            diff_waits_problem_preprobe_unpredict, correlation_plot_cell_unpredict{3}, 'Prediction accuracy and amount of learning')
        axis([-25 25 -3 8]); ylabel('Last day - first day (d`)')
        var_name = 'PredictAccuracy_AmountLearn_corr'; 
    	%print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
    plot2var(diff_waits_problem_preprobe_predict, correlation_plot_cell_predict{4}, ...
            diff_waits_problem_preprobe_unpredict, correlation_plot_cell_unpredict{4}, 'Prediction accuracy and memory stability')
        axis([-25 25 -.6 .8]); ylabel('Correlation between pre and post probe (r)')
        var_name = 'PredictAccuracy_MemoryStability_corr'; 
    	%print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
    plot2var(diff_waits_problem_preprobe_predict, correlation_plot_cell_predict{5}, ...
            diff_waits_problem_preprobe_unpredict, correlation_plot_cell_unpredict{5}, 'Prediction accuracy and final discrimination')
        axis([-25 25 -2 9]); ylabel('Correlation between pre and post probe (r)')
        var_name = 'PredictAccuracy_Lastday_corr'; 
    	print([figure_path var_name], '-dpdf', '-painters', '-bestfit')

%mixed model prep
num_subjects = 80; num_problems = 6;
subj_num_mtx = repmat((1:num_subjects)', 1, num_problems); % subject matrix
stage_num_mtx = repmat(1:6, num_subjects, 1); % stage matrix
model_str = 'var1~var2+(1|subject)';


% mixed models
lme_first_day_discrim = lme_function(model_str, first_day_discrimination, accuracy_of_prediction, subj_num_mtx)
lme_days2crit = lme_function(model_str, days_to_criterion, accuracy_of_prediction, subj_num_mtx)
lme_amount = lme_function(model_str, amount_learning, accuracy_of_prediction, subj_num_mtx)
lme_stability = lme_function(model_str, probe_over_probe_corr, accuracy_of_prediction, subj_num_mtx)
lme_last_day_discrim = lme_function(model_str, last_day_discrimination, accuracy_of_prediction, subj_num_mtx)

function lme = lme_function(model_str, var1, var2, subj_num_mtx)
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
    tbl = table(var1(:), var2(:), subj_num_mtx(:),'VariableNames',{inputname(2), inputname(3), 'subject'});
    tbl.subject = categorical(tbl.subject);
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
    hold on; plot(xlim, [0 0], 'k--')
    hold on; plot([0 0], ylim, 'k--')

end
    
    
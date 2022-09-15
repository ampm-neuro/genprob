green_color = [21 154 72]./255; blue_color =  [43 73 159]./255;
light_green_color = [162 201 129]./255; light_blue_color = [102 164 213]./255;
figure_path = 'C:\Users\ampm1\Documents\manuscripts\gentrain\manuscript_figures\sfig01\';


%% probe 0, problem 1, probe 1

figure; 

% probe 0
subplot(1,3,1); hold on    
    
    % unpredict probe curve
    probe_paths = get_file_paths_targeted('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_hivar', 'preprobe_01');
    plot_allprobes(probe_paths); 

    % predict probe curve
    probe_paths = get_file_paths_targeted('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar', 'preprobe_01');
    plot_allprobes(probe_paths);

    % aesthetics
    rich_bounds_prob('mevar', 1);
    ylim([0 25]); ylabel('Wait durations (s)')
    xlim([4500 38500]); xticks([5000 35000]); xticklabels({'5', '35'}); xlabel('Freq. (kHz)')
    axis square
    title('Probe 0')

% problem 1
subplot(1,3,2); hold on

    % predict
    plot_all_subj_problem('train_mevar', 1, 1, 2)
    plot_all_subj_problem('train_mevar', 1, 2, 2)

    % unpredict
    plot_all_subj_problem('train_hivar', 1, 1, 2)
    plot_all_subj_problem('train_hivar', 1, 2, 2)
    
    % aesthetics
    title('Problem 1')
    ylim([0 25]); ylabel('Wait durations (s)')
    xlim([0 42]); xticks([1 41]); xticklabels({'5', '35'}); xlabel('Freq. (kHz)')
    rich_bounds_prob('mevar', 1, 1);
    
% probe 1
subplot(1,3,3); hold on    
    
    % unpredict probe curve
    probe_paths = get_file_paths_targeted('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_hivar', 'preprobe_02');
    plot_allprobes(probe_paths); 

    % predict probe curve
    probe_paths = get_file_paths_targeted('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar', 'preprobe_02');
    plot_allprobes(probe_paths);

    % aesthetics
    rich_bounds_prob('mevar', 1);
    ylim([0 25]); ylabel('Wait durations (s)')
    xlim([4500 38500]); xticks([5000 35000]); xticklabels({'5', '35'}); xlabel('Freq. (kHz)')
    axis square
    title('Probe 1')

% save
set(gcf, 'Position', [34 134 1648 497])
var_name = 'Probe0Problem1Probe1'; 
%print(['C:\Users\ampm1\Documents\manuscripts\gentrain\manuscript_figures\Sfig01\' var_name], '-dpdf', '-painters', '-bestfit')



%% probe 5, problem 6, probe 6
figure; 

% probe 5
subplot(1,3,1); hold on    
    
    % unpredict probe curve
    probe_paths = get_file_paths_targeted('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_hivar', 'preprobe_06');
    plot_allprobes(probe_paths); 

    % predict probe curve
    probe_paths = get_file_paths_targeted('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar', 'preprobe_06');
    plot_allprobes(probe_paths);

    % aesthetics
    rich_bounds_prob('mevar', 6);
    ylim([0 25]); ylabel('Wait durations (s)')
    xlim([4500 38500]); xticks([5000 35000]); xticklabels({'5', '35'}); xlabel('Freq. (kHz)')
    axis square
    title('Probe 5')

% problem 1
subplot(1,3,2); hold on

    % predict
    plot_all_subj_problem('train_mevar', 6, 1, 2)
    plot_all_subj_problem('train_mevar', 6, 2, 2)

    % unpredict
    plot_all_subj_problem('train_hivar', 6, 1, 2)
    plot_all_subj_problem('train_hivar', 6, 2, 2)
    
    % aesthetics
    title('Problem 6')
    ylim([0 25]); ylabel('Wait durations (s)')
    xlim([0 42]); xticks([1 41]); xticklabels({'5', '35'}); xlabel('Freq. (kHz)')
    rich_bounds_prob('mevar', 6, 1);
    
% probe 1
subplot(1,3,3); hold on    
    
    % unpredict probe curve
    probe_paths = get_file_paths_targeted('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_hivar', 'postprobe_01');
    plot_allprobes(probe_paths); 

    % predict probe curve
    probe_paths = get_file_paths_targeted('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar', 'postprobe_01');
    plot_allprobes(probe_paths);

    % aesthetics
    rich_bounds_prob('mevar', 6);
    ylim([0 25]); ylabel('Wait durations (s)')
    xlim([4500 38500]); xticks([5000 35000]); xticklabels({'5', '35'}); xlabel('Freq. (kHz)')
    axis square
    title('Probe 6')

% save
set(gcf, 'Position', [34 134 1648 497])
var_name = 'Probe5Problem6Probe6'; 
%print(['C:\Users\ampm1\Documents\manuscripts\gentrain\manuscript_figures\Sfig01\' var_name], '-dpdf', '-painters', '-bestfit')



%% rewarded trial wait times histogram
window_size = 2;
bins = 2:window_size:60;
data_path_mevar = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar';
probe_sessions_only_mevar = get_file_paths_targeted(data_path_mevar, 'probe');
data_path_hivar = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_hivar';
probe_sessions_only_hivar = get_file_paths_targeted(data_path_hivar, 'probe');
combine_probe_sessions = [probe_sessions_only_mevar;probe_sessions_only_hivar];
combine_probe_sessions = combine_probe_sessions(~contains(combine_probe_sessions, 'notone'));
[all_trl_mtx, session_number, problem_number] = ALL_trl_mtx(combine_probe_sessions);
rwd_ct = nan(1, length(bins)-1); abandon_ct = nan(size(rwd_ct)); rwd_prop = nan(size(rwd_ct)); 
for itb = 1:length(bins)-1
    delay_duration_idx = all_trl_mtx(:,4)>bins(itb) & all_trl_mtx(:,4)<=bins(itb+1);
    rwd_ct(itb) = sum(all_trl_mtx(all_trl_mtx(:,3)==1 & ~isnan(all_trl_mtx(:,11)) & delay_duration_idx));
    abandon_ct(itb) = sum(all_trl_mtx(all_trl_mtx(:,3)==1 & isnan(all_trl_mtx(:,11)) & delay_duration_idx));
    rwd_prop(itb) = rwd_ct(itb) / (rwd_ct(itb) + abandon_ct(itb));
end
stack_plot_input = [rwd_ct; abandon_ct]'; stack_plot_input = stack_plot_input./sum(stack_plot_input(:));
figure; hold on; 
bar(bins(1:end-1)+(window_size/2), stack_plot_input, 'stacked')
set(gca,'TickLength',[0, 0]);
axis square
xlim([0 60]); ylim([0 .26])
plot([2 2], ylim, 'k--')
ylabel('Probability')
xlabel('Wait time (s)')
title('Reward available trial wait times (probes only)')
legend({'Rewarded', 'Abandoned'})

    % inset
    figure; hold on;
    plot(bins(1:end-1)+(window_size/2), rwd_prop)
    plot(bins(1:end-1)+(window_size/2), 1-rwd_prop)
    set(gca,'TickLength',[0, 0]);
    axis square
    xlabel('Delay duration (s)')
    ylabel('Proportion')
    plot([2 2], ylim, 'k--')
    title('Rewards received or abandoned')


    
%% nonrewarded trial wait times histogram
no_reward_trial_waits = all_trl_mtx(all_trl_mtx(:,12)<60 & all_trl_mtx(:,12)>0 & all_trl_mtx(:,3)==0, 12)+2;
figure; hold on; histogram(no_reward_trial_waits, 0:2:60, 'normalization', 'probability')
median_norwd = median(no_reward_trial_waits);
set(gca,'TickLength',[0, 0]);
axis square
xlim([0 60]); ylim([0 .16])
plot([2 2], ylim, 'k--'); ylabel('Probability')
xlabel('Wait time (s)')
title('Reward unnavailable trial wait times (probes only)')

    % inset
    figure; hold on;
    plot(bins(1:end-1)+(window_size/2), cumsum(histcounts(no_reward_trial_waits, bins(bins>=2)))./length(no_reward_trial_waits))
    set(gca,'TickLength',[0, 0]);
    axis square
    xlabel('Delay duration (s)')
    ylabel('Cumulative proportion')
    plot([2 2], ylim, 'k--')
    title('Reward unnavailable trial wait times')

    
    
%% example problem session
%{
% search for examples
data_path_mevar = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar';
subject_folders = get_folder_paths_all(data_path_mevar);
% iterate through subjects
for isubj = 1:size(subject_folders,1)
    % load last session
    problem6 = get_file_paths_targeted(subject_folders{isubj}, 'gen12');
    load(problem6{end}, 'trl_mtx')
    % plot if 100 trials
    if size(trl_mtx,1)>=100
        try
        figure; plot_trials_trlmtx(trl_mtx); title(num2str(isubj)); drawnow
        catch; disp('catch'); end
    end
end
%}

% trials plot
load('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar\677761m4\gen11_mevar05-03.mat', 'trl_mtx')
figure; plot_trials_trlmtx(trl_mtx); title('677761m4, problem5 last');
set(gca, 'YDir','reverse'); yticks([1 10:10:100]); axis([-3 60 0.5 101]);
xlabel('Time from head entry (s)')

% no reward rich vs poor
figure; wait_times_plot(trl_mtx,3); ylim([0 60]); xlim([4500 38500])



%% example probe sessions
% search for examples
%{
data_path_mevar = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_hivar';
subject_folders = get_folder_paths_all(data_path_mevar);
% iterate through subjects
for isubj = 1:size(subject_folders,1)
    % load last session
    probe6 = get_file_paths_targeted(subject_folders{isubj}, 'preprobe_02');
    load(probe6{1}, 'trl_mtx')
    % plot if 100 trials
    if size(trl_mtx,1)>=82
        try
            figure; 
            subplot(1,2,1)
            plot_trials_trlmtx(trl_mtx); 
            xlim([-5 60])
            set(gca, 'YDir','reverse')
            yticks([1 10:10:70 82])
            subplot(1,2,2); hold on
            wait_times_plot(trl_mtx); wait_times_plot_tonecolor(trl_mtx);
            ylim([0 60])
        catch; disp('catch'); end
        title(num2str(isubj));
    end
    drawnow
end
%}

% trials plot
%load('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar\689831m3\gen14_postprobe_01_01d.mat', 'trl_mtx')
%figure; 
%plot_trials_trlmtx(trl_mtx); title('689831m3, postprobe01');
load('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_hivar\691359m4\gen14_preprobe_02_01d.mat', 'trl_mtx')
figure;
plot_trials_trlmtx(trl_mtx); title('691359m4, preprobe02');
xlim([-3 60])
ylim([0.5 83])
set(gca, 'YDir','reverse')
yticks([1 10:10:70 82])

% curve plot
figure; hold on
wait_times_plot(trl_mtx); wait_times_plot_tonecolor(trl_mtx); 
ylim([0 60]); xlim([4500 38500])


%% compare how task execution changed over learning

data_path_mevar = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar';
data_path_hivar = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_hivar';
%data_path_mevar = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_imaging_hpc';
%data_path_hivar = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_hivar_imaging_hpc';


% preallocate
np_times_predict = nan(40,7);
he_times_predict = nan(40,7);
np_times_unpredict = nan(40,7);
he_times_unpredict = nan(40,7);
trial_ct_predict = nan(40,7);
trial_ct_unpredict = nan(40,7);
session_duration_predict = nan(40,7);
session_duration_unpredict = nan(40,7);

% iterate through probe sessions
for iprobe = 1:7
    
    % relevant probe sessions
    if iprobe < 7 
        probe_sessions_only_mevar = get_file_paths_targeted(data_path_mevar, ['preprobe_0' num2str(iprobe)]);
        probe_sessions_only_hivar = get_file_paths_targeted(data_path_hivar, ['preprobe_0' num2str(iprobe)]);
    else
        probe_sessions_only_mevar = get_file_paths_targeted(data_path_mevar, 'postprobe_01');
        probe_sessions_only_hivar = get_file_paths_targeted(data_path_hivar, 'postprobe_01');
    end
    
    % predict
        
        % iterate through sessions
        for isubj = 1:size(probe_sessions_only_mevar, 1)
            
            % load
            load(probe_sessions_only_mevar{isubj}, 'trl_mtx')
            
            % nose poke duration
            np_durations = trl_mtx(:,9)-trl_mtx(:,6);
            np_durations(np_durations<0 | np_durations>10) = nan;
            np_times_predict(isubj, iprobe) = nanmean(np_durations);
    
            % time to head entry
            tt_he = trl_mtx(:,10)-trl_mtx(:,9);
            tt_he(tt_he<0 | tt_he>4) = nan;
            he_times_predict(isubj, iprobe) = nanmean(tt_he);
            
            % number of trials
            trial_ct_predict(isubj, iprobe) = size(trl_mtx,1);
            
            % session duration
            session_duration_predict(isubj, iprobe) = trl_mtx(end,1)./60;
            
        end
        
        
    % unpredict
        
        % iterate through sessions
        for isubj = 1:size(probe_sessions_only_hivar, 1)
            
            % load
            load(probe_sessions_only_hivar{isubj}, 'trl_mtx')
            
            % nose poke duration
            np_durations = trl_mtx(:,9)-trl_mtx(:,6);
            np_durations(np_durations<0 | np_durations>10) = nan;
            np_times_unpredict(isubj, iprobe) = nanmean(np_durations);
    
            % time to head entry
            tt_he = trl_mtx(:,10)-trl_mtx(:,9);
            tt_he(tt_he<0 | tt_he>4) = nan;
            he_times_unpredict(isubj, iprobe) = nanmean(tt_he);
            
            % number of trials
            trial_ct_unpredict(isubj, iprobe) = size(trl_mtx,1);
            
            % session duration
            session_duration_unpredict(isubj, iprobe) = trl_mtx(end,1)./60;
            
        end
end

% plot
%
% nose poke
%{
figure; hold on;
errorbar_mtx(cat(3, np_times_predict, np_times_unpredict), [blue_color; green_color], [light_blue_color; light_green_color]);
title('Nose poke durations')
ylabel('Time (s)')
xticks(1:7); xticklabels(0:6); xlabel('Probe')
axis square; set(gca,'TickLength',[0, 0]);
var_name = 'Nose_poke_duration';
print(['C:\Users\ampm1\Documents\manuscripts\gentrain\manuscript_figures\Sfig01\' var_name], '-dpdf', '-painters', '-bestfit')


        %  stats: mixed model test for interaction
        % repeated measures anova
        datamtx = [np_times_unpredict; np_times_predict];
        between_subjs = [zeros(40,1); ones(40,1)];
        rmtbl = simple_mixed_anova(datamtx, between_subjs, {'problem'}, {'group'});
        ranova_interact = [rmtbl.DF(3) rmtbl.DF(5) rmtbl.F(5) rmtbl.pValue(5)];
%}

% head entry initiation time
%{
figure; hold on
errorbar_mtx(cat(3, he_times_unpredict, he_times_predict), [blue_color; green_color], [light_blue_color; light_green_color]);
title('Time to head entry')
ylabel('Time (seconds)')
xticks(1:7); xticklabels(0:6); xlabel('Probe')
axis square; set(gca,'TickLength',[0, 0]);
var_name = 'NP_to_HE_time';
print(['C:\Users\ampm1\Documents\manuscripts\gentrain\manuscript_figures\Sfig01\' var_name], '-dpdf', '-painters', '-bestfit')

        %  stats: mixed model test for interaction
        % repeated measures anova
        datamtx = [he_times_unpredict; he_times_predict];
        between_subjs = [zeros(40,1); ones(40,1)];
        rmtbl = simple_mixed_anova(datamtx, between_subjs, {'problem'}, {'group'});
        ranova_interact = [rmtbl.DF(3) rmtbl.DF(5) rmtbl.F(5) rmtbl.pValue(5)];
%}        
        
% nose poke + head entry = trial initiation time
figure; hold on
errorbar_mtx(cat(3, np_times_unpredict+he_times_unpredict, np_times_predict+he_times_predict), [blue_color; green_color], [light_blue_color; light_green_color]);
title('Trial initiation')
ylabel('Time (seconds)')
xticks(1:7); xticklabels(0:6); xlabel('Probe')
axis square; set(gca,'TickLength',[0, 0]); ylim([0 4])
var_name = 'Trial_initiation_time';
print(['C:\Users\ampm1\Documents\manuscripts\gentrain\manuscript_figures\Sfig01\' var_name], '-dpdf', '-painters', '-bestfit')

        %  stats: mixed model test for interaction
        % repeated measures anova
        datamtx = [np_times_unpredict+he_times_unpredict; np_times_predict+he_times_predict];
        between_subjs = [zeros(40,1); ones(40,1)];
        rmtbl = simple_mixed_anova(datamtx, between_subjs, {'problem'}, {'group'});
        ranova_interact = [rmtbl.DF(3) rmtbl.DF(5) rmtbl.F(5) rmtbl.pValue(5)];


% session duration
figure; hold on
errorbar_mtx(cat(3, session_duration_unpredict, session_duration_predict), [blue_color; green_color], [light_blue_color; light_green_color]);
title('Session duration')
ylabel('Time (minutes)'); ylim([0 61])
xticks(1:7); xticklabels(0:6); xlabel('Probe')
axis square; set(gca,'TickLength',[0, 0]);
var_name = 'Session_duration';
print(['C:\Users\ampm1\Documents\manuscripts\gentrain\manuscript_figures\Sfig01\' var_name], '-dpdf', '-painters', '-bestfit')

        %  stats: mixed model test for interaction
        % repeated measures anova
        datamtx = [session_duration_unpredict; session_duration_predict];
        between_subjs = [zeros(40,1); ones(40,1)];
        rmtbl = simple_mixed_anova(datamtx, between_subjs, {'problem'}, {'group'});
        ranova_interact = [rmtbl.DF(3) rmtbl.DF(5) rmtbl.F(5) rmtbl.pValue(5)];
        
%
figure; hold on;
errorbar_mtx(cat(3, trial_ct_unpredict, trial_ct_predict), [blue_color; green_color], [light_blue_color; light_green_color])
title('Trial count')
ylabel('Time (s)'); ylim([0 83])
xticks(1:7); xticklabels(0:6); xlabel('Probe')
axis square; set(gca,'TickLength',[0, 0]);
hold on; plot(xlim, [1 1].*41, 'k--')
var_name = 'Trial_count';
print(['C:\Users\ampm1\Documents\manuscripts\gentrain\manuscript_figures\Sfig01\' var_name], '-dpdf', '-painters', '-bestfit')
%}

        %  stats: mixed model test for interaction
        % repeated measures anova
        % THIS ISNT NORMALLY DISTRIBUTED OBVI
        datamtx = [trial_ct_unpredict; trial_ct_predict];
        between_subjs = [zeros(40,1); ones(40,1)];
        rmtbl = simple_mixed_anova(datamtx, between_subjs, {'problem'}, {'group'});
        ranova_interact = [rmtbl.DF(3) rmtbl.DF(5) rmtbl.F(5) rmtbl.pValue(5)];

        
        
%% Compare probe behavior correlation matrices (see figure 2)

[probe_waits_predict, ~, cm_17_predict, offdiag_cell_predict, comp2one_cell_predict, comp2seven_cell_predict] = probe_corr_mtx('train_mevar');
close
[probe_waits_unpredict, ~, cm_17_unpredict, offdiag_cell_unpredict, comp2one_cell_unpredict, comp2seven_cell_unpredict] = probe_corr_mtx('train_hivar');
close

% compare similarity of responses on all probes to probe 6 
% (test degree to which nearer-in-time probes are similar to each other)
%
figure; hold on; 
errorbar_plot(comp2seven_cell_unpredict, 1, [], light_blue_color, blue_color);
errorbar_plot(comp2seven_cell_predict, 1, [], light_green_color, green_color);
xlabel('Probe')
xticks(1:6)
xticklabels(0:5);
ylabel('Similarity (r)')
plot(xlim, [1 1].*0, 'k--')
ylim([-1 1])
title('Compare to probe 7')
axis square
var_name = 'ProbeBeh_comp2seven_comp'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')

        %  stats: repeated measures anova for interaction
        datamtx = [cell2mat(comp2seven_cell_unpredict); cell2mat(comp2seven_cell_predict)];
        between_subjs = [zeros(40,1); ones(40,1)];
        rmtbl = simple_mixed_anova(datamtx, between_subjs, {'problem'}, {'group'});
        ranova_comp2seven_stats_comp = [rmtbl.DF(3) rmtbl.DF(5) rmtbl.F(5) rmtbl.pValue(5)];
    
        % between group ttests
        pvals_comp2seven_comp = nan(1,6); stats_cell_comp2seven_comp = cell(1,6);
        for itest = 1:6
            [~, pvals_comp2seven_comp(itest), ~, stats_cell_comp2seven_comp{itest}] = ttest2(comp2seven_cell_unpredict{itest}, comp2seven_cell_predict{itest});
        end

%}

% compare similarity of responses on all probes to probe 0 
% (test degree to which nearer-in-time probes are similar to each other)
% main effect of group (predict is slightly below 0), no interaction
%{
figure; hold on; 
errorbar_plot(comp2one_cell_unpredict, 1, [], light_blue_color, blue_color);
errorbar_plot(comp2one_cell_predict, 1, [], light_green_color, green_color);
plot(xlim, [1 1].*0, 'k--')
xticks(1:6)
xlabel('Probe')
ylabel('Similarity (r)')
ylim([-1 1])
title('Compare to probe 1')
axis square
var_name = 'ProbeBeh_comp2one_comp'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')

        %  stats: repeated measures anova for interaction
        datamtx = [cell2mat(comp2one_cell_unpredict); cell2mat(comp2one_cell_predict)];
        between_subjs = [zeros(40,1); ones(40,1)];
        rmtbl = simple_mixed_anova(datamtx, between_subjs, {'problem'}, {'group'});
        ranova_comp2one_stats_comp = [rmtbl.DF(3) rmtbl.DF(5) rmtbl.F(5) rmtbl.pValue(5)];
    
        % between group ttests
        pvals_comp2one_comp = nan(1,6); stats_cell_comp2one_comp = cell(1,6);
        for itest = 1:6
            [~, pvals_comp2one_comp(itest), ~, stats_cell_comp2one_comp{itest}] = ttest2(comp2one_cell_unpredict{itest}, comp2one_cell_predict{itest});
        end
        
%}
        
        
%% Fit curves and evaluate normal and logistic components over training
% see compare_model_responses script

% fit curves to every probe
%[logist_coef_cell_corrected_mevar, norm_coef_cell_mevar] = ALL_probe_vs_behavior_fit('train_mevar');
%[logist_coef_cell_corrected_hivar, norm_coef_cell_hivar] = ALL_probe_vs_behavior_fit('train_hivar');
%save('coefs', 'logist_coef_cell_corrected_mevar', 'norm_coef_cell_mevar', 'logist_coef_cell_corrected_hivar', 'norm_coef_cell_hivar')
load('coefs', 'logist_coef_cell_corrected_mevar', 'norm_coef_cell_mevar', 'logist_coef_cell_corrected_hivar', 'norm_coef_cell_hivar')

% probes of interest
poi = [1 2 6 7];

% logistic curves
figure; hold on
errorbar_plot(logist_coef_cell_corrected_hivar(poi), 1, [], light_blue_color, blue_color)
errorbar_plot(logist_coef_cell_corrected_mevar(poi), 1, [], light_green_color, green_color)
plot(xlim, [1 1].*0, 'k--')
xticks(1:4); xticklabels([0 1 5 6])
ylim([-20 30]); ylabel('Coefficient (s)')
xlabel('Probe')
title('Logistic')
var_name = 'Logistic_curve';
print(['C:\Users\ampm1\Documents\manuscripts\gentrain\manuscript_figures\Sfig01\' var_name], '-dpdf', '-painters', '-bestfit')

        %  stats: mixed model test for interaction
        % repeated measures anova
        datamtx = [cell2mat(logist_coef_cell_corrected_hivar(poi)); cell2mat(logist_coef_cell_corrected_mevar(poi))];
        between_subjs = [zeros(40,1); ones(40,1)];
        rmtbl = simple_mixed_anova(datamtx, between_subjs, {'problem'}, {'group'});
        ranova_interact = [rmtbl.DF(3) rmtbl.DF(5) rmtbl.F(5) rmtbl.pValue(5)];
        
        % ttest between groups
        ttest_pvals_logistic = nan(1,length(poi));
        ttest_stats_logistic = cell(1,length(poi));
        for ipoi = 1:length(poi)
            [~,ttest_pvals_logistic(ipoi),~,ttest_stats_logistic{ipoi}] = ttest2(cell2mat(logist_coef_cell_corrected_hivar(poi(ipoi))),cell2mat(logist_coef_cell_corrected_mevar(poi(ipoi))));
        end
        

% normal curves
figure; hold on
errorbar_plot(norm_coef_cell_hivar(poi), 1, [], light_blue_color, blue_color)
errorbar_plot(norm_coef_cell_mevar(poi), 1, [], light_green_color, green_color)
plot(xlim, [1 1].*0, 'k--')
xticks(1:4); xticklabels([0 1 5 6])
ylim([-20 30]); ylabel('Coefficient (s)')
xlabel('Probe')
title('Normal')
var_name = 'Normal_curve';
print(['C:\Users\ampm1\Documents\manuscripts\gentrain\manuscript_figures\Sfig01\' var_name], '-dpdf', '-painters', '-bestfit')

        %  stats: mixed model test for interaction
        % repeated measures anova
        datamtx = [cell2mat(norm_coef_cell_hivar(poi)); cell2mat(norm_coef_cell_mevar(poi))];
        between_subjs = [zeros(40,1); ones(40,1)];
        rmtbl = simple_mixed_anova(datamtx, between_subjs, {'problem'}, {'group'});
        ranova_interact = [rmtbl.DF(3) rmtbl.DF(5) rmtbl.F(5) rmtbl.pValue(5)];

        % ttest between groups
        ttest_pvals_normal = nan(1,length(poi));
        ttest_stats_normal = cell(1,length(poi));
        for ipoi = 1:length(poi)
            [~,ttest_pvals_normal(ipoi),~,ttest_stats_normal{ipoi}] = ttest2(cell2mat(norm_coef_cell_hivar(poi(ipoi))),cell2mat(norm_coef_cell_mevar(poi(ipoi))));
        end


%% context test
tbl = table(colon_op([wait_times_all{1}; wait_times_all{2}]), colon_op(repmat(1:41, 16,1)), colon_op([zeros(8,41);ones(8,41)]), colon_op(repmat([1:8 1:8]', 1,41)), 'VariableNames', {'Waits', 'Tones', 'Context', 'Subjects'});
tbl.Context = categorical(tbl.Context);
tbl.Subjects = categorical(tbl.Subjects);

[p,tbl,stats,terms] = anovan(tbl.Waits,{tbl.Tones,tbl.Context,tbl.Subjects},'model',3,'random',3,'varnames',{'Tones', 'Context', 'Subjects'});






















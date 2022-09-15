% a directory of code used to create figures for labmeet6
% be careful to select the right subjects
% this code varies between using all subejcts and only subjects that
% completed every problem and probe

%% training problem order
rwd_tones_2t_var2

%% single subject figures

% trial by trial plots
plot_trials_trlmtx_altcolor(trl_mtx) % colored by rich and poor tones
plot_trials_trlmtx(trl_mtx) % colored by tone frequency

% full learning curve (dropped extra/second sessions)
continuous_learning(11, '664773m1');
continuous_learning_singlesesh(session_path, 11); % see also

% single probe sessions
figure; wait_times_plot(trl_mtx,3);
hold on; [coefEsts] = plot_normal_fit_subj(trl_mtx, 3);

% probe and problem overlays
probe_problem_interact_plots('670035m3', 'train_hivar')

% subject summaries
ALL_subj_sum('train_hivar');
subj_summary('660469m2')


%% average learning curve from all subjects normalized
[all_prefs, subj_idx, problem_idx, sesh_idx] = ALL_mean_cont_learn_proportionate('train_hivar');
[all_prefs, subj_idx, problem_idx, sesh_idx] = ALL_mean_cont_learn('train_hivar'); % see also


%%  early late preference for rich tone
% BE SURE TO SELECT WHAT PROPORTION OF LEARNING CURVE COUNTS AS EARLY/LATE
[all_prefs, subj_idx, problem_idx, first_mtx, last_mtx, all_sesh_ct, ...
    learn_rate_mtx, delta_mtx, min_discrim_mtx, unq_subjects_beh, ...
    first_mtx_rich, first_mtx_poor, last_mtx_rich, last_mtx_poor] = ALL_cont_learn_delta_dot(1, 'train_hivar');


%% rich tone preference during probes
folder_in = 'mevar_fin';
[rich_waits, unq_subjects_probe_rich] = probe_wait_times(folder_in, 1:7, 16:26);


%% days to criterion
df = 'two_tone\train_hivar'; 
[all_stage_learn_mtx, all_stage_learn_cell] = ALL_stage_learn(df, 1:6, 2, 1);


%% correlations

% early discrim VS trials to criterion
% early discrim VS late discrim
% BE SURE TO SELECT WHAT PROPORTION OF LEARNING CURVE COUNTS AS EARLY/LATE
[all_prefs, subj_idx, problem_idx, first_mtx, last_mtx, all_sesh_ct, ...
    learn_rate_mtx, delta_mtx, min_discrim_mtx, unq_subjects_beh, ...
    first_mtx_rich, first_mtx_poor, last_mtx_rich, last_mtx_poor] = ALL_cont_learn_delta_dot(1, 'train_hivar');

% probe behavior VS training problem behavior
ALL_probe_vs_behavior('train_hivar') % ALSO HAS MIXED MODELS
% see also: ALL_probe_vs_behavior_fit

% poor tone problem change VS poor tone probe change
% rich tone problem change VS rich tone probe change
ALL_probe_vs_behavior('train_hivar')

% just first and last day (or, you know, whatever day)
doi_MevarHivar


%% all probe average responses
for i = 1:8
    figure;
    plot_probe_stages(i, 'mevar'); 
    ylim([0 25])
end

%% plotting problem behavior, one dot per mouse
plot_problemlearning_bysubj


%% preference for rich tone freqs during probe, line
compare_type = 2;
figure; hold on
[rich_waits, unq_subjects_probe_rich] = probe_wait_times('train_hivar', 1:8, 16:26, compare_type);
[rich_waits, unq_subjects_probe_rich] = probe_wait_times('train_mevar', 1:8, 16:26, compare_type);
hold on; plot(xlim, [1 1].*0, 'k--')
title('mevar rich zone')

compare_type = 2;
figure; hold on
[rich_waits, unq_subjects_probe_rich] = probe_wait_times('train_hivar', 1:8, [16 18 21 24 26], compare_type);
[rich_waits, unq_subjects_probe_rich] = probe_wait_times('train_mevar', 1:8, [16 18 21 24 26], compare_type);
hold on; plot(xlim, [1 1].*0, 'k--')
title('mevar rich tones')

compare_type = 2;
figure; hold on
[rich_waits, unq_subjects_probe_rich] = probe_wait_times('train_hivar', 1:8, [5 13 21 29 37], compare_type);
[rich_waits, unq_subjects_probe_rich] = probe_wait_times('train_mevar', 1:8, [5 13 21 29 37], compare_type);
hold on; plot(xlim, [1 1].*0, 'k--')
title('hivar rich tones')



%% normalized probe over probe change in wait times
ALL_probe_vs_behavior('train_hivar')

%% behavior correlation matrix
df = 'two_tone\train_mevar_fin'; [bcm, bcc] = beh_corr_mtx(df);

%% probe correlation matrix (GREAT!)
probe_corr_mtx('train_mevar')
probe_corr_mtx('train_hivar')



%% plot probes

% compare hi and med variance
for i = 1:8
    figure; hold on; 
    plot_probe_stages(i, 'train_hivar_fin'); 
    rich_bounds_prob('hivar', i-1);
    plot_probe_stages(i, 'train_mevar_fin');
    rich_bounds_prob('mevar', i-1);
    ylim([0 30])
    if i < 8
        title(['After problem ' num2str(i-1)]); 
    else
        title(['One week later']); 
        rich_bounds_prob('hivar', i-1);
    end
    ylim([0 40])
    axis square
end

% overlay all of a single var type
figure; hold on; for i = 1:7; plot_probe_stages(i, 'hivar'); end

%% Probe specificity
% how do wait times to trained tones compare to nearby ones during probe
probe_specificity('mevar_fin', 2, [], 2);

%% compare logistic, normal coefs between conditions
compare_model_responses


%% probe over probe changes
ALL_probe_over_probe_delta_comb




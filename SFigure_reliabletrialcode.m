green_color = [21 154 72]./255; blue_color =  [43 73 159]./255;
light_green_color = [162 201 129]./255; light_blue_color = [102 164 213]./255;
figure_path = 'C:\Users\ampm1\Documents\manuscripts\gentrain\manuscript_figures\sfig05\';
predict_subjs = [{'651049m1'} {'658648m2'} {'691359m1'} {'690330m1'} {'691359m2'} {'696944m3'} {'699437m4'} {'699438m3'}];
unpredict_subjs = [{'687034m4'} {'683472m2'} {'699437m2'} {'699438m1'} {'683472m3'}];

% for time warping
tses = [0.2 1.0 1.0 2.0 3.0 2.0];

%% example time warping
load('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_imaging_hpc\691359m2\LED_gen12_mevar06-02.mat')
neuron_num = 7;
figure; hold on
tw_activity_plot_trial_hm(trl_mtx, trl_idx, medass_cell, frame_times, traces, neuron_num, 1:size(trl_mtx,1), 1, [5 25]);
xlim([490 3000]); axis square
var_name = 'time_warp_eg_raw'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')

figure; hold on
tw_activity_plot_trial_hm_warp(trl_mtx, trl_idx, frame_times, traces, neuron_num, 1:size(trl_mtx,1), ceil([sum(tses(1:3)) sum(tses(4:end))]), tses);
xlim([75 Inf]); set(gca, 'YDir','normal'); axis square
var_name = 'time_warp_eg_warp'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')


%% reliable sequential population activity from trial-to-trial
%{
[within_session_corr_means_predict, session_cell_means_predict] = within_session_reliability_frfields(predict_subjs, 1:1:19, tses);
[within_session_corr_means_unpredict, session_cell_means_unpredict] = within_session_reliability_frfields(unpredict_subjs, 1:1:19, tses);
save('within_session_trialcorrs_all', 'within_session_corr_means_predict', 'within_session_corr_means_unpredict', 'session_cell_means_predict', 'session_cell_means_unpredict')
%}
load('within_session_trialcorrs_all', 'within_session_corr_means_predict', 'within_session_corr_means_unpredict', 'session_cell_means_predict', 'session_cell_means_unpredict')
within_session_corr_means_predict_probesonly = within_session_corr_means_predict(:,1:3:19);
within_session_corr_means_unpredict_probesonly = within_session_corr_means_unpredict(:,1:3:19);

% plot
figure; hold on
errorbar_mtx(within_session_corr_means_unpredict_probesonly, blue_color, light_blue_color)
errorbar_mtx(within_session_corr_means_predict_probesonly, green_color, light_green_color)
xticks(1:7); xticklabels(0:6); xlabel('Probe')
ylabel('Mean trial-by-trial activity correlation (r)')
axis square
var_name = 'within_session_trialcorrs'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')

    datamtx_WithinStability_probes = [within_session_corr_means_unpredict_probesonly; within_session_corr_means_predict_probesonly];
    num_subjects = length(predict_subjs) + length(unpredict_subjs); 
    num_sessions = size(within_session_corr_means_unpredict_probesonly,2);
    subj_num_mtx = repmat((1:num_subjects)', 1, num_sessions); % subject matrix
    session_num_mtx = repmat(1:num_sessions, num_subjects, 1); % stage matrix
    group_mtx = [zeros(length(predict_subjs), num_sessions); ones(length(unpredict_subjs), num_sessions)];
    model_str = 'PopOverlap~Session*Group+(1|Subject)';
    tbl = table(datamtx_WithinStability_probes(:), session_num_mtx(:), group_mtx(:), subj_num_mtx(:), 'VariableNames', {'PopOverlap', 'Session', 'Group', 'Subject'});
    tbl.Subject = categorical(tbl.Subject);
    lme_PreprobeProblem = fitlme(tbl, model_str)


    
%% psuedo session example showing reliable firing from trial to trial

% all probe sessions
%{
[combine_trial_activity_mtx_predict, subject_ids_predict, session_nums_predict, cell_regist_mtx_cell_predict] = image_trial_activity_mtx_multisubj(predict_subjs, 1:3:19, 1:41, tses);
[combine_trial_activity_mtx_unpredict, subject_ids_unpredict, session_nums_unpredict, cell_regist_mtx_cell_unpredict] = image_trial_activity_mtx_multisubj(unpredict_subjs, 1:3:19, 1:41, tses);
save('combine_trial_activity_predict', 'combine_trial_activity_mtx_predict', 'subject_ids_predict', 'session_nums_predict', 'cell_regist_mtx_cell_predict')
save('combine_trial_activity_unpredict', 'combine_trial_activity_mtx_unpredict', 'subject_ids_unpredict', 'session_nums_unpredict', 'cell_regist_mtx_cell_unpredict')
%}
load('combine_trial_activity_predict'); load('combine_trial_activity_unpredict')

% get sort indices
[~,sort_idx_predict] = sort_rows_by_peak(norm_mtx(nanmean(combine_trial_activity_mtx_predict,3)')');
[~,sort_idx_unpredict] = sort_rows_by_peak(norm_mtx(nanmean(combine_trial_activity_mtx_unpredict,3)')');

% plot trials
toi = 1:5:41;
figure; 
for itrl = 1:length(toi)
    subplot(1, length(toi), itrl); hold on
    imagesc(norm_mtx(combine_trial_activity_mtx_predict(sort_idx_predict,:,toi(itrl))')')
    set(gca, 'YDir','reverse'); ylim([0.5 size(combine_trial_activity_mtx_predict,1)+0.5])
    title(['Trial ' num2str(toi(itrl))])
    yticklabels([])
    for ief = event_frame(1:end-1) 
       plot((ief-min(event_frame)+1).*[1 1], ylim, 'r-', 'linewidth', 1.0) 
    end
    xlim([0.5 event_frame(4)-event_frame(1)+1.5])
end
sgtitle('Predictable training')
set(gcf, 'Position', [478 245 1245 298])
var_name = 'within_session_reliable_fig_predict'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')

figure; 
for itrl = 1:length(toi)
    subplot(1, length(toi), itrl); hold on
    imagesc(norm_mtx(combine_trial_activity_mtx_unpredict(sort_idx_unpredict,:,toi(itrl))')')
    set(gca, 'YDir','reverse'); ylim([0.5 size(combine_trial_activity_mtx_unpredict,1)+0.5])
    title(['Trial ' num2str(toi(itrl))])
    yticklabels([])
    for ief = event_frame 
       plot((ief-min(event_frame)+1).*[1 1], ylim, 'r-', 'linewidth', 1.0) 
    end
    xlim([0.5 event_frame(4)-event_frame(1)+1.5])
end
sgtitle('Unpredictable training')
set(gcf, 'Position', [478 245 1245 298])
var_name = 'within_session_reliable_fig_unpredict'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')

    
%% low dimensional population trajectory plots


%% specificity changes with reactivation


%% reliability changes with reactivation







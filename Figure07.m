% do reactivation scores correlate with information stored in final probe
% response


%% information in every subject's final probe about each training problem

% get waits from every probe to every problem
%{
[diff_waits_problem_predict, ratio_waits_problem_predict] = ALL_probe_vs_behavior_terpolate('train_mevar_imaging_hpc');
[diff_waits_problem_unpredict, ratio_waits_problem_unpredict] = ALL_probe_vs_behavior_terpolate('train_hivar_imaging_hpc');
save('dif_waits_to_all_probs_nofit', 'diff_waits_problem_predict', 'diff_waits_problem_unpredict', 'ratio_waits_problem_predict', 'ratio_waits_problem_unpredict')
%}
load('dif_waits_to_all_probs_nofit', 'diff_waits_problem_predict', 'diff_waits_problem_unpredict', 'ratio_waits_problem_predict', 'ratio_waits_problem_unpredict')

% just final probes
diff_waits_problem_predict = squeeze(diff_waits_problem_predict(7,:,:))';
diff_waits_problem_unpredict = squeeze(diff_waits_problem_unpredict(7,:,:))';
ratio_waits_problem_predict = squeeze(ratio_waits_problem_predict(7,:,:))';
ratio_waits_problem_unpredict = squeeze(ratio_waits_problem_unpredict(7,:,:))';

% reactivation score during each 
load('rscores_raw', 'rscore_predict_PreprobeProblem', 'rscore_unpredict_PreprobeProblem', 'rscore_predict_ProblemFinalprobe',...
    'rscore_unpredict_ProblemFinalprobe', 'rscore_predict_PreprobeFinalprobe', 'rscore_unpredict_PreprobeFinalprobe')
%load('rscores_timecorrect', 'rscore_predict_PreprobeProblem', 'rscore_unpredict_PreprobeProblem', 'rscore_predict_ProblemFinalprobe',...
%    'rscore_unpredict_ProblemFinalprobe', 'rscore_predict_PreprobeFinalprobe', 'rscore_unpredict_PreprobeFinalprobe')


% plot
figure
subplot(1,3,1); hold on
plot2var(rscore_predict_PreprobeProblem, diff_waits_problem_predict, rscore_unpredict_PreprobeProblem, diff_waits_problem_unpredict, 'Preprobe reactivation during problem')
ylabel('Rich wait - poor wait (s)'); xlabel('Raw rscore')
subplot(1,3,2); hold on
plot2var(rscore_predict_ProblemFinalprobe, diff_waits_problem_predict, rscore_unpredict_ProblemFinalprobe, diff_waits_problem_unpredict, 'Problem reactivation during final probe')
ylabel('Rich wait - poor wait (s)'); xlabel('Raw rscore')
subplot(1,3,3); hold on
plot2var(rscore_predict_PreprobeFinalprobe, diff_waits_problem_predict, rscore_unpredict_PreprobeFinalprobe, diff_waits_problem_unpredict, 'Preprobe reactivation during final probe')
ylabel('Rich wait - poor wait (s)'); xlabel('Raw rscore')
sgtitle('Diff waits')

figure
subplot(1,3,1); hold on
plot2var(rscore_predict_PreprobeProblem, ratio_waits_problem_predict, rscore_unpredict_PreprobeProblem, ratio_waits_problem_unpredict, 'Preprobe reactivation during problem')
ylabel('Rich-poor / rich+poor (s)'); xlabel('Raw rscore')
subplot(1,3,2); hold on
plot2var(rscore_predict_ProblemFinalprobe, ratio_waits_problem_predict, rscore_unpredict_ProblemFinalprobe, ratio_waits_problem_unpredict, 'Problem reactivation during final probe')
ylabel('Rich-poor / rich+poor (s)'); xlabel('Raw rscore')
subplot(1,3,3); hold on
plot2var(rscore_predict_PreprobeFinalprobe, ratio_waits_problem_predict, rscore_unpredict_PreprobeFinalprobe, ratio_waits_problem_unpredict, 'Preprobe reactivation during final probe')
ylabel('Rich-poor / rich+poor (s)'); xlabel('Raw rscore')
sgtitle('Ratio waits')


%% model

datamtx_adj_problem = [rscore_predict_PreprobeProblem; rscore_unpredict_PreprobeProblem];
datamtx_last_problem = [rscore_predict_ProblemFinalprobe; rscore_unpredict_ProblemFinalprobe];
datamtx_last_probes = [rscore_predict_PreprobeFinalprobe; rscore_unpredict_PreprobeFinalprobe];
diff_waits_all = [diff_waits_problem_predict; diff_waits_problem_unpredict];
ratio_waits_all = [ratio_waits_problem_predict; ratio_waits_problem_unpredict];
num_subjects = size(rscore_predict_PreprobeProblem,1) + size(rscore_unpredict_PreprobeProblem,1); 
num_sessions = size(rscore_predict_PreprobeProblem,2);
subj_num_mtx = repmat((1:num_subjects)', 1, num_sessions); % subject matrix
session_num_mtx = repmat(1:num_sessions, num_subjects, 1); % stage matrix
group_mtx = [zeros(size(rscore_predict_PreprobeProblem,1), num_sessions); ones(size(rscore_unpredict_PreprobeProblem,1), num_sessions)];

model_str = 'Waits~Rscore+(1|Subject)';
tbl_diff_TrainingReact = table(diff_waits_all(:), rscore_adj_problem(:), subj_num_mtx(:), 'VariableNames', {'Waits', 'Rscore', 'Subject'}); tbl.Subject = categorical(tbl.Subject);
tbl_ratio_TrainingReact = table(ratio_waits_all(:), rscore_adj_problem(:), subj_num_mtx(:), 'VariableNames', {'Waits', 'Rscore', 'Subject'}); tbl.Subject = categorical(tbl.Subject);
tbl_diff_ProblemreactFin = table(diff_waits_all(:), rscore_last_problem(:), subj_num_mtx(:), 'VariableNames', {'Waits', 'Rscore', 'Subject'}); tbl.Subject = categorical(tbl.Subject);
tbl_ratio_ProblemreactFin = table(ratio_waits_all(:), rscore_last_problem(:), subj_num_mtx(:), 'VariableNames', {'Waits', 'Rscore', 'Subject'}); tbl.Subject = categorical(tbl.Subject);
tbl_diff_PreprobereactFin = table(diff_waits_all(:), rscore_last_probes(:), subj_num_mtx(:), 'VariableNames', {'Waits', 'Rscore', 'Subject'}); tbl.Subject = categorical(tbl.Subject);
tbl_ratio_PreprobereactFin = table(ratio_waits_all(:), rscore_last_probes(:), subj_num_mtx(:), 'VariableNames', {'Waits', 'Rscore', 'Subject'}); tbl.Subject = categorical(tbl.Subject);

lme_PreprobeProblem_diff = fitlme(tbl_diff_TrainingReact, model_str)
lme_PreprobeProblem_ratio = fitlme(tbl_ratio_TrainingReact, model_str)
lme_ProblemFinprobe_diff = fitlme(tbl_diff_TrainingReact, model_str)
lme_ProblemFinprobe_ratio = fitlme(tbl_ratio_TrainingReact, model_str)
lme_PreprobeFinprobe_diff = fitlme(tbl_diff_TrainingReact, model_str)
lme_PreprobeFinprobe_ratio = fitlme(tbl_ratio_TrainingReact, model_str)




%% internal functions
function plot2var(var1_predict, var2_predict, var1_unpredict, var2_unpredict, title_str)
% vars need to be matrices with (subjects, samples)

    % colormap
    green_colors = [162 221 154; 118 205 117; 065 181 093; 037 149 069; 007 120 056; 013 078 032]./255;
    blue_colors = [158 202 225; 106 173 214; 067 146 198; 035 114 181; 009 083 157; 033 052 104]./255;
    
    % open figure
    %figure; hold on
    
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
    axis square
end
green_color = [21 154 72]./255; blue_color =  [43 73 159]./255;
light_green_color = [162 201 129]./255; light_blue_color = [102 164 213]./255;
figure_path = 'C:\Users\ampm1\Documents\manuscripts\gentrain\manuscript_figures\fig06\';
predict_subjs = [{'651049m1'} {'658648m2'} {'691359m1'} {'690330m1'} {'691359m2'} {'696944m3'} {'699437m4'} {'699438m3'}];
unpredict_subjs = [{'687034m4'} {'683472m2'} {'699437m2'} {'699438m1'} {'683472m3'}];

% for time warping
tses = [0.2 1.0 1.0 2.0 3.0 2.0];
event_frame = cumsum(tses(1:end-1)).*100; %from second np to reward delivery/ quit

% key imaging sessions
preprobe_session_numbers = 1:3:16;
firstday_session_numbers = 2:3:17;
lastday_session_numbers = 3:3:18;
%% behavior
%{
[~, ~, diff_waits_problem_preprobe_predict, diff_waits_problem_postprobe_predict, correlation_plot_cell_predict] = ALL_probe_vs_behavior_fit('train_mevar_imaging_hpc');
[~, ~, diff_waits_problem_preprobe_unpredict, diff_waits_problem_postprobe_unpredict, correlation_plot_cell_unpredict] = ALL_probe_vs_behavior_fit('train_hivar_imaging_hpc');
save('sfigure02_2', 'diff_waits_problem_preprobe_predict', 'diff_waits_problem_preprobe_unpredict', 'correlation_plot_cell_predict', 'correlation_plot_cell_unpredict')
%}
load('sfigure02_2', 'diff_waits_problem_preprobe_predict', 'diff_waits_problem_preprobe_unpredict', 'correlation_plot_cell_predict', 'correlation_plot_cell_unpredict')



%% Reinstatement score predicted by accuracy of retrieval

% overlap reinstatement during learing
load('overlap_raw', 'prop_predict', 'prop_unpredict')

% activity pattern reinstatement during learning
load('Popcorr_raw', 'subj_corr_means_mtx_predict', 'subj_corr_means_mtx_unpredict')

% reinstatement score
rscore_predict = prop_predict + (prop_predict.*subj_corr_means_mtx_predict);
rscore_unpredict = prop_unpredict + (prop_unpredict.*subj_corr_means_mtx_unpredict);

% remove identity
for isesh = 1:size(rscore_predict,1)
    rscore_predict(isesh, isesh, :) = nan;
    rscore_unpredict(isesh, isesh, :) = nan;
end



    
    
%% model comparing overlap matrices
subj_corr_means_mtx_predict_noduplicates = rscore_predict;
subj_corr_means_mtx_unpredict_noduplicates = rscore_unpredict;
for isesh1 = 1:size(rscore_predict,1)
    subj_corr_means_mtx_predict_noduplicates(isesh1, isesh1, :) = nan;
    subj_corr_means_mtx_unpredict_noduplicates(isesh1, isesh1, :) = nan;
end
datamtx = [subj_corr_means_mtx_predict_noduplicates(:); subj_corr_means_mtx_unpredict_noduplicates(:)];

num_subjects = length(predict_subjs)+length(unpredict_subjs); 
num_sessions1 = 19;
num_sessions2 = 19;
subj_num_mtx = repmat((1:num_subjects), num_sessions1*num_sessions2, 1); subj_num_mtx = subj_num_mtx(:); % subject matrix
session_num_mtx1 = repmat((1:num_sessions1)', num_sessions1, num_subjects); session_num_mtx1 = session_num_mtx1(:); % stage matrix
session_num_mtx2 = repmat(1:num_sessions2, num_sessions2, num_subjects); session_num_mtx2 = session_num_mtx2(:); % stage matrix
group_mtx = [zeros(length(predict_subjs)*num_sessions1*num_sessions2, 1); ones(length(unpredict_subjs)*num_sessions1*num_sessions2, 1)];

nnan_idx = ~isnan(datamtx);
datamtx = datamtx(nnan_idx);
subj_num_mtx = subj_num_mtx(nnan_idx);
session_num_mtx1 = session_num_mtx1(nnan_idx);
session_num_mtx2 = session_num_mtx2(nnan_idx);
group_mtx = group_mtx(nnan_idx);

model_str = 'MemorySimilarity~Session1*Session2*Group+(1|Subject)';
tbl = table(datamtx, session_num_mtx1, session_num_mtx2, group_mtx, subj_num_mtx, 'VariableNames', {'MemorySimilarity', 'Session1', 'Session2', 'Group', 'Subject'});
tbl.Subject = categorical(tbl.Subject);
lme_ReactivationMatrices = fitlme(tbl, model_str)


    
%% compute values of interest
    
preprobe_session_numbers = 1:3:16;
firstday_session_numbers = 2:3:17;
lastday_session_numbers = 3:3:18;
postprobe_session_numbers = 4:3:19;

% preallcate subjects x probes
rscore_predict_PreprobeProblem = nan(length(predict_subjs), length(preprobe_session_numbers));
rscore_predict_ProblemFinalprobe = nan(length(predict_subjs), length(preprobe_session_numbers));
rscore_predict_PreprobeFinalprobe = nan(length(predict_subjs), length(preprobe_session_numbers));
rscore_unpredict_PreprobeProblem = nan(length(unpredict_subjs), length(preprobe_session_numbers));
rscore_unpredict_ProblemFinalprobe = nan(length(unpredict_subjs), length(preprobe_session_numbers));
rscore_unpredict_PreprobeFinalprobe = nan(length(unpredict_subjs), length(preprobe_session_numbers));

% iterate through sessions
for isesh = 1:length(preprobe_session_numbers)
    
    % adjacent probe sessions
    %
        % iterate through predict subjects
        for isubj_predict = 1:length(predict_subjs)
            rscore_predict_PreprobeProblem(isubj_predict, isesh) = ...
                mean(rscore_predict(preprobe_session_numbers(isesh), [firstday_session_numbers(isesh) lastday_session_numbers(isesh)], isubj_predict));
        end
        % iterate through unpredict subjects
        for isubj_unpredict = 1:length(unpredict_subjs)
            rscore_unpredict_PreprobeProblem(isubj_unpredict, isesh) = ...
                mean(rscore_unpredict(preprobe_session_numbers(isesh), [firstday_session_numbers(isesh) lastday_session_numbers(isesh)], isubj_unpredict));
        end
        
    % compare all problems with last probe
    %
        % iterate through predict subjects
        for isubj_predict = 1:length(predict_subjs)
            rscore_predict_ProblemFinalprobe(isubj_predict, isesh) = ...
                mean([rscore_predict(firstday_session_numbers(isesh), end, isubj_predict);... 
                      rscore_predict(lastday_session_numbers(isesh), end, isubj_predict)]); % compared to probe 6
        end
        % iterate through unpredict subjects
        for isubj_unpredict = 1:length(unpredict_subjs)
            rscore_unpredict_ProblemFinalprobe(isubj_unpredict, isesh) = ...
                mean([rscore_unpredict(firstday_session_numbers(isesh), end, isubj_unpredict);...
                      rscore_unpredict(lastday_session_numbers(isesh), end, isubj_unpredict)]);% compared to probe 6
        end
        
    % compare all probes with last probe
    %
        % iterate through predict subjects
        for isubj_predict = 1:length(predict_subjs)
            rscore_predict_PreprobeFinalprobe(isubj_predict, isesh) = rscore_predict(preprobe_session_numbers(isesh), postprobe_session_numbers(end), isubj_predict);
        end
        % iterate through unpredict subjects
        for isubj_unpredict = 1:length(unpredict_subjs)
            rscore_unpredict_PreprobeFinalprobe(isubj_unpredict, isesh) = rscore_unpredict(preprobe_session_numbers(isesh), postprobe_session_numbers(end), isubj_unpredict);
        end
end

%% PLOTS

% combine data for later mixed models
datamtx_PreprobeProblem = [rscore_predict_PreprobeProblem; rscore_unpredict_PreprobeProblem];
    datamtx_PreprobeProblem_common = datamtx_PreprobeProblem(:, [1 end]);
datamtx_PreprobeFinalprobe = [rscore_predict_PreprobeFinalprobe; rscore_unpredict_PreprobeFinalprobe];
    datamtx_PreprobeFinalprobe_common = datamtx_PreprobeFinalprobe(:, [1 end]);
datamtx_ProblemFinalprobe = [rscore_predict_ProblemFinalprobe; rscore_unpredict_ProblemFinalprobe];
    datamtx_ProblemFinalprobe_common = datamtx_ProblemFinalprobe(:, [1 end]);

% plot preprobe vs problem
%
f1 = figure;
% predict
subplot(1,4,1); hold on;
errorbar_plot(num2cell(rscore_predict_PreprobeProblem,1), 1, [], light_green_color, green_color)
xticks(1:size(rscore_predict_PreprobeProblem,2)); xlabel('Probe number'); xticklabels(0:size(rscore_predict_PreprobeFinalprobe,2)-1);
ylim([-0.1 1.1]); yticks(0:.2:1); ylabel('Reactivation score')
axis square; title('Predict')
% unpredict
subplot(1,4,2); hold on;
errorbar_plot(num2cell(rscore_unpredict_PreprobeProblem,1), 1, [], light_blue_color, blue_color)
xticks(1:size(rscore_predict_PreprobeProblem,2)); xlabel('Probe number'); xticklabels(0:size(rscore_predict_PreprobeFinalprobe,2)-1);
ylim([-0.1 1.1]); yticks(0:.2:1); ylabel('Reactivation score')
axis square; title('Unpredict')
% direct compare first/last
subplot(1,4,3); hold on;
errorbar_plot(num2cell(rscore_unpredict_PreprobeProblem(:,[1 6]),1), 1, [], light_blue_color, blue_color)
errorbar_plot(num2cell(rscore_predict_PreprobeProblem(:,[1 6]),1), 1, [], light_green_color, green_color)
xticks(1:size(rscore_predict_PreprobeProblem,2)); xlabel('Probe number'); xticklabels([0 5]);
ylim([-0.1 1.1]); yticks(0:.2:1); ylabel('Reactivation score')
axis square; title('Compare')
sgtitle('Preprobe and problem')

    % mixed model full
    num_subjects = length(predict_subjs) + length(unpredict_subjs); 
    num_sessions = length(preprobe_session_numbers);
    subj_num_mtx = repmat((1:num_subjects)', 1, num_sessions); % subject matrix
    session_num_mtx = repmat(1:num_sessions, num_subjects, 1); % stage matrix
    group_mtx = [zeros(length(predict_subjs), num_sessions); ones(length(unpredict_subjs), num_sessions)];
    model_str = 'PopOverlap~Session*Group+(1|Subject)';
    tbl = table(datamtx_PreprobeProblem(:), session_num_mtx(:), group_mtx(:), subj_num_mtx(:), 'VariableNames', {'PopOverlap', 'Session', 'Group', 'Subject'});
    tbl.Subject = categorical(tbl.Subject);
    lme_PreprobeProblem_errorbar = fitlme(tbl, model_str)

        % one way anova
        % fill in nan with group condition mean
        rscore_predict_PreprobeProblem_anovahold = rscore_predict_PreprobeProblem;
        for icond = 1:size(rscore_predict_PreprobeProblem_anovahold,2)
            rscore_predict_PreprobeProblem_anovahold(isnan(rscore_predict_PreprobeProblem_anovahold(:,icond)),icond) = nanmean(rscore_predict_PreprobeProblem_anovahold(:,icond));
        end
        rmtbl_adj_PreprobeProblem_predict = simple_mixed_anova(rscore_predict_PreprobeProblem_anovahold);
        rmtbl_adj_PreprobeProblem_unpredict = simple_mixed_anova(rscore_unpredict_PreprobeProblem);
        
    % ttest compare to zero
    pval_predict_PreprobeProblem_0 = nan(size(rscore_predict_PreprobeProblem,2),1); 
        stats_predict_PreprobeProblem_0 = cell(size(rscore_predict_PreprobeProblem,2),1);
    pval_unpredict_PreprobeProblem_0 = nan(size(rscore_unpredict_PreprobeProblem,2),1); 
        stats_unpredict_PreprobeProblem_0 = cell(size(rscore_unpredict_PreprobeProblem,2),1);
    for icomp = 1:size(rscore_predict_PreprobeProblem,2)
        [~,pval_predict_PreprobeProblem_0(icomp),~,stats_predict_PreprobeProblem_0{icomp}] = ttest(rscore_predict_PreprobeProblem(:,icomp));
        [~,pval_unpredict_PreprobeProblem_0(icomp),~,stats_unpredict_PreprobeProblem_0{icomp}] = ttest(rscore_unpredict_PreprobeProblem(:,icomp));
    end

    % ttest stats within groups
    pval_adj_predict_firstlast = nan(size(rscore_predict_PreprobeProblem,2),1); stats_adj_predict_firstlast = cell(size(rscore_predict_PreprobeProblem,2),1);
    pval_adj_unpredict_firstlast = nan(size(rscore_unpredict_PreprobeProblem,2),1); stats_adj_unpredict_firstlast = cell(size(rscore_unpredict_PreprobeProblem,2),1);
    for icomp = 2:size(rscore_predict_PreprobeProblem,2)
        [~,pval_adj_predict_firstlast(icomp),~,stats_adj_predict_firstlast{icomp}] = ttest(rscore_predict_PreprobeProblem(:,1), rscore_predict_PreprobeProblem(:,icomp));
        [~,pval_adj_unpredict_firstlast(icomp),~,stats_adj_unpredict_firstlast{icomp}] = ttest(rscore_unpredict_PreprobeProblem(:,1), rscore_unpredict_PreprobeProblem(:,icomp));
    end

    % ttest stats between groups
    pval_adj_comp_firstlast = nan(size(rscore_unpredict_PreprobeProblem,2),1); stats_adj_comppval_adj_comp_firstlast = cell(size(rscore_unpredict_PreprobeProblem,2),1);
    for icomp = 1:size(rscore_predict_PreprobeProblem,2)
        [~,pval_adj_comp_firstlast(icomp),~,stats_adj_comppval_adj_comp_firstlast{icomp}] = ttest2(rscore_predict_PreprobeProblem(:,icomp), rscore_unpredict_PreprobeProblem(:,icomp));
    end
    xpos = 0;
    for icomp = [1 size(rscore_predict_PreprobeProblem,2)]
        xpos = xpos + 1;
        pval_text(pval_adj_comp_firstlast(icomp), xpos-.4, 1)
    end
    
    % mixed model common
    num_subjects = length(predict_subjs) + length(unpredict_subjs); 
    num_sessions = 2;
    subj_num_mtx = repmat((1:num_subjects)', 1, num_sessions); % subject matrix
    session_num_mtx = repmat(1:num_sessions, num_subjects, 1); % stage matrix
    group_mtx = [zeros(length(predict_subjs), num_sessions); ones(length(unpredict_subjs), num_sessions)];
    model_str = 'PopOverlap~Session*Group+(1|Subject)';
    tbl = table(datamtx_PreprobeProblem_common(:), session_num_mtx(:), group_mtx(:), subj_num_mtx(:), 'VariableNames', {'PopOverlap', 'Session', 'Group', 'Subject'});
    tbl.Subject = categorical(tbl.Subject);
    lme_PreprobeProblem_common = fitlme(tbl, model_str)

% plot problem vs probe 6
%
f3 = figure;
% predict
subplot(1,4,1); hold on;
errorbar_plot(num2cell(rscore_predict_ProblemFinalprobe,1), 1, [], light_green_color, green_color)
xticks(1:size(rscore_predict_ProblemFinalprobe,2)); xlabel('Problem number'); xticklabels(1:6);
ylim([-0.1 1.1]); yticks(0:.2:1); ylabel('Reactivation score')
axis square; title('Predictable')
% unpredict
subplot(1,4,2); hold on;
errorbar_plot(num2cell(rscore_unpredict_ProblemFinalprobe,1), 1, [], light_blue_color, blue_color)
xticks(1:size(rscore_unpredict_ProblemFinalprobe,2)); xlabel('Problem number'); xticklabels(1:6);
ylim([-0.1 1.1]); yticks(0:.2:1); ylabel('Reactivation score')
axis square; title('Unpredictable')
% direct compare
subplot(1,4,3); hold on;
errorbar_plot(num2cell(rscore_unpredict_ProblemFinalprobe(:,[1 6]),1), 1, [], light_blue_color, blue_color)
errorbar_plot(num2cell(rscore_predict_ProblemFinalprobe(:,[1 6]),1), 1, [], light_green_color, green_color)
xticks(1:size(rscore_predict_ProblemFinalprobe,2)); xlabel('Problem number'); xticklabels([1 6]);
ylim([-0.1 1.1]); yticks(0:.2:1); ylabel('Reactivation score')
axis square; title('Direct compare')
sgtitle('problem and final probe')
 
    % mixed model
    num_subjects = length(predict_subjs) + length(unpredict_subjs); 
    num_sessions = length(preprobe_session_numbers);
    subj_num_mtx = repmat((1:num_subjects)', 1, num_sessions); % subject matrix
    session_num_mtx = repmat(1:num_sessions, num_subjects, 1); % stage matrix
    group_mtx = [zeros(length(predict_subjs), num_sessions); ones(length(unpredict_subjs), num_sessions)];
    model_str = 'PopOverlap~Session*Group+(1|Subject)';
    tbl = table(datamtx_ProblemFinalprobe(:), session_num_mtx(:), group_mtx(:), subj_num_mtx(:), 'VariableNames', {'PopOverlap', 'Session', 'Group', 'Subject'});
    tbl.Subject = categorical(tbl.Subject);
    lme_ProblemFinalprobe = fitlme(tbl, model_str)

        % one way anova
        % fill in nan with group condition mean
        rscore_predict_ProblemFinalprobe_anovahold = rscore_predict_ProblemFinalprobe;
        for icond = 1:size(rscore_predict_ProblemFinalprobe_anovahold,2)
            rscore_predict_ProblemFinalprobe_anovahold(isnan(rscore_predict_ProblemFinalprobe_anovahold(:,icond)),icond) = nanmean(rscore_predict_ProblemFinalprobe_anovahold(:,icond));
        end
        rmtbl_adj_ProblemFinalprobe_predict = simple_mixed_anova(rscore_predict_ProblemFinalprobe_anovahold);
        rmtbl_adj_ProblemFinalprobe_unpredict = simple_mixed_anova(rscore_unpredict_ProblemFinalprobe);

    % ttest compare to zero
    pval_last_predict_firstlast_0 = nan(size(rscore_predict_ProblemFinalprobe,2),1); 
        stats_last_predict_firstlast_0 = cell(size(rscore_predict_ProblemFinalprobe,2),1);
    pval_last_unpredict_firstlast_0 = nan(size(rscore_unpredict_ProblemFinalprobe,2),1); 
        stats_unpredict_firstlast_0 = cell(size(rscore_unpredict_ProblemFinalprobe,2),1);
    for icomp = 1:size(rscore_predict_ProblemFinalprobe,2)
        [~,pval_last_predict_firstlast_0(icomp),~,stats_last_predict_firstlast_0{icomp}] = ttest(rscore_predict_ProblemFinalprobe(:,icomp));
        [~,pval_last_unpredict_firstlast_0(icomp),~,stats_unpredict_firstlast_0{icomp}] = ttest(rscore_unpredict_ProblemFinalprobe(:,icomp));
    end

    % ttest stats within groups
    pval_last_predict_firstlast = nan(size(rscore_predict_ProblemFinalprobe,2),1); 
        stats_last_predict_firstlast = cell(size(rscore_predict_ProblemFinalprobe,2),1);
    pval_last_unpredict_firstlast = nan(size(rscore_unpredict_ProblemFinalprobe,2),1); 
        stats_last_unpredict_firstlast = cell(size(rscore_unpredict_ProblemFinalprobe,2),1);
    for icomp = 2:size(rscore_predict_ProblemFinalprobe,2)
        [~,pval_last_predict_firstlast(icomp),~,stats_last_predict_firstlast{icomp}] = ttest(rscore_predict_ProblemFinalprobe(:,1), rscore_predict_ProblemFinalprobe(:,icomp));
        [~,pval_last_unpredict_firstlast(icomp),~,stats_last_unpredict_firstlast{icomp}] = ttest(rscore_unpredict_ProblemFinalprobe(:,1), rscore_unpredict_ProblemFinalprobe(:,icomp));
    end

    % ttest stats between group
    pval_last_comp_firstlast = nan(size(rscore_unpredict_ProblemFinalprobe,2),1); stats_last_comp_firstlast = cell(size(rscore_unpredict_ProblemFinalprobe,2),1);
    for icomp = 1:size(rscore_predict_ProblemFinalprobe,2)
        [~,pval_last_comp_firstlast(icomp),~,stats_last_comp_firstlast{icomp}] = ttest2(rscore_predict_ProblemFinalprobe(:,icomp), rscore_unpredict_ProblemFinalprobe(:,icomp));
    end
    xpos = 0;
    for icomp = [1 size(rscore_unpredict_ProblemFinalprobe,2)]
        xpos = xpos + 1;
        pval_text(pval_last_comp_firstlast(icomp), xpos-.4, 1)
    end
    
    % mixed model common
    num_subjects = length(predict_subjs) + length(unpredict_subjs); 
    num_sessions = 2;
    subj_num_mtx = repmat((1:num_subjects)', 1, num_sessions); % subject matrix
    session_num_mtx = repmat(1:num_sessions, num_subjects, 1); % stage matrix
    group_mtx = [zeros(length(predict_subjs), num_sessions); ones(length(unpredict_subjs), num_sessions)];
    model_str = 'PopOverlap~Session*Group+(1|Subject)';
    tbl = table(datamtx_ProblemFinalprobe_common(:), session_num_mtx(:), group_mtx(:), subj_num_mtx(:), 'VariableNames', {'PopOverlap', 'Session', 'Group', 'Subject'});
    tbl.Subject = categorical(tbl.Subject);
    lme_ProblemFinalprobe_common = fitlme(tbl, model_str)
    

% plot preprobe vs probe 6
%
f2 = figure;
% predict
subplot(1,4,1); hold on;
errorbar_plot(num2cell(rscore_predict_PreprobeFinalprobe,1), 1, [], light_green_color, green_color)
xticks(1:size(rscore_predict_PreprobeProblem,2)); xlabel('Probe number'); xticklabels(0:size(rscore_predict_PreprobeFinalprobe,2)-1);
ylim([-0.1 1.1]); yticks(0:.2:1); ylabel('Reactivation score')
axis square; title('Predict')
% unpredict
subplot(1,4,2); hold on;
errorbar_plot(num2cell(rscore_unpredict_PreprobeFinalprobe,1), 1, [], light_blue_color, blue_color)
xticks(1:size(rscore_predict_PreprobeProblem,2)); xlabel('Probe number'); xticklabels(0:size(rscore_predict_PreprobeFinalprobe,2)-1);
ylim([-0.1 1.1]); yticks(0:.2:1); ylabel('Reactivation score')
axis square; title('Unpredict')
% direct compare
subplot(1,4,3); hold on;
errorbar_plot(num2cell(rscore_unpredict_PreprobeFinalprobe(:,[1 6]),1), 1, [], light_blue_color, blue_color)
errorbar_plot(num2cell(rscore_predict_PreprobeFinalprobe(:,[1 6]),1), 1, [], light_green_color, green_color)
xticks(1:size(rscore_predict_PreprobeProblem,2)); xlabel('Probe number'); xticklabels([0 5]);
ylim([-0.1 1.1]); yticks(0:.2:1); ylabel('Reactivation score')
axis square; title('Compare')
sgtitle('Preprobe and final probe')



    % mixed model 
    num_subjects = length(predict_subjs) + length(unpredict_subjs); 
    num_sessions = length(preprobe_session_numbers);
    subj_num_mtx = repmat((1:num_subjects)', 1, num_sessions); % subject matrix
    session_num_mtx = repmat(1:num_sessions, num_subjects, 1); % stage matrix
    group_mtx = [zeros(length(predict_subjs), num_sessions); ones(length(unpredict_subjs), num_sessions)];
    model_str = 'PopOverlap~Session*Group+(1|Subject)';
    tbl = table(datamtx_PreprobeFinalprobe(:), session_num_mtx(:), group_mtx(:), subj_num_mtx(:), 'VariableNames', {'PopOverlap', 'Session', 'Group', 'Subject'});
    tbl.Subject = categorical(tbl.Subject);
    lme_PreprobeFinalprobe = fitlme(tbl, model_str)
    
        % one way anova
        % fill in nan with group condition mean
        rscore_predict_PreprobeFinalprobe_anovahold = rscore_predict_PreprobeFinalprobe;
        for icond = 1:size(rscore_predict_PreprobeFinalprobe_anovahold,2)
            rscore_predict_PreprobeFinalprobe_anovahold(isnan(rscore_predict_PreprobeFinalprobe_anovahold(:,icond)),icond) = nanmean(rscore_predict_PreprobeFinalprobe_anovahold(:,icond));
        end
        rmtbl_adj_PreprobeFinalprobe_predict = simple_mixed_anova(rscore_predict_PreprobeFinalprobe_anovahold);
        rmtbl_adj_PreprobeFinalprobe_unpredict = simple_mixed_anova(rscore_unpredict_PreprobeFinalprobe);

    % ttest compare to zero
    pval_last_predict_0 = nan(size(rscore_predict_PreprobeFinalprobe,2),1); stats_last_predict_0 = cell(size(rscore_predict_PreprobeFinalprobe,2),1);
    pval_last_unpredict_0 = nan(size(rscore_unpredict_PreprobeFinalprobe,2),1); stats_last_unpredict_0 = cell(size(rscore_unpredict_PreprobeFinalprobe,2),1);
    for icomp = 1:size(rscore_predict_PreprobeFinalprobe,2)
        [~,pval_last_predict_0(icomp),~,stats_last_predict_0{icomp}] = ttest(rscore_predict_PreprobeFinalprobe(:,icomp));
        [~,pval_last_unpredict_0(icomp),~,stats_last_unpredict_0{icomp}] = ttest(rscore_unpredict_PreprobeFinalprobe(:,icomp));
    end

    % ttest stats within groups
    pval_last_predict = nan(size(rscore_predict_PreprobeFinalprobe,2),1); stats_last_predict = cell(size(rscore_predict_PreprobeFinalprobe,2),1);
    pval_last_unpredict = nan(size(rscore_unpredict_PreprobeFinalprobe,2),1); stats_last_unpredict = cell(size(rscore_unpredict_PreprobeFinalprobe,2),1);
    for icomp = 2:size(rscore_predict_PreprobeFinalprobe,2)
        [~,pval_last_predict(icomp),~,stats_last_predict{icomp}] = ttest(rscore_predict_PreprobeFinalprobe(:,1), rscore_predict_PreprobeFinalprobe(:,icomp));
        [~,pval_last_unpredict(icomp),~,stats_last_unpredict{icomp}] = ttest(rscore_unpredict_PreprobeFinalprobe(:,1), rscore_unpredict_PreprobeFinalprobe(:,icomp));
    end

    % ttest stats between group
    pval_last_comp = nan(size(rscore_unpredict_PreprobeFinalprobe,2),1); stats_last_comp = cell(size(rscore_unpredict_PreprobeFinalprobe,2),1);
    for icomp = 1:size(rscore_predict_PreprobeFinalprobe,2)
        [~,pval_last_comp(icomp),~,stats_last_comp{icomp}] = ttest2(rscore_predict_PreprobeFinalprobe(:,icomp), rscore_unpredict_PreprobeFinalprobe(:,icomp));
    end
    xpos = 0;
    for icomp = [1 size(rscore_predict_PreprobeFinalprobe,2)]
        xpos = xpos + 1;
        pval_text(pval_last_comp(icomp), xpos-.4, 1)
    end
    
    % mixed model common
    num_subjects = length(predict_subjs) + length(unpredict_subjs); 
    num_sessions = 2;
    subj_num_mtx = repmat((1:num_subjects)', 1, num_sessions); % subject matrix
    session_num_mtx = repmat(1:num_sessions, num_subjects, 1); % stage matrix
    group_mtx = [zeros(length(predict_subjs), num_sessions); ones(length(unpredict_subjs), num_sessions)];
    model_str = 'PopOverlap~Session*Group+(1|Subject)';
    tbl = table(datamtx_PreprobeFinalprobe_common(:), session_num_mtx(:), group_mtx(:), subj_num_mtx(:), 'VariableNames', {'PopOverlap', 'Session', 'Group', 'Subject'});
    tbl.Subject = categorical(tbl.Subject);
    lme_PreprobeFinalprobe_common = fitlme(tbl, model_str)

    
save('rscores_raw', 'rscore_predict_PreprobeProblem', 'rscore_unpredict_PreprobeProblem', 'rscore_predict_ProblemFinalprobe',...
    'rscore_unpredict_ProblemFinalprobe', 'rscore_predict_PreprobeFinalprobe', 'rscore_unpredict_PreprobeFinalprobe')


    
    
% Overlap comparisons (correlations)
%
%[~, ~, diff_waits_problem_preprobe_predict, diff_waits_problem_postprobe_predict, correlation_plot_cell_predict] = ALL_probe_vs_behavior_fit('train_mevar_imaging_hpc');
%for ifig = 1:34; close; end
%[~, ~, diff_waits_problem_preprobe_unpredict, diff_waits_problem_postprobe_unpredict, correlation_plot_cell_unpredict] = ALL_probe_vs_behavior_fit('train_hivar_imaging_hpc');
%for ifig = 1:34; close; end
load('sfigure02_2', 'diff_waits_problem_preprobe_predict', 'diff_waits_problem_preprobe_unpredict', 'correlation_plot_cell_predict', 'correlation_plot_cell_unpredict')

accuracy_of_prediction = [diff_waits_problem_preprobe_predict; diff_waits_problem_preprobe_unpredict];
first_day_discrimination = [correlation_plot_cell_predict{1}; correlation_plot_cell_unpredict{1}];
% correlation_plot_cell{1} = first day discrim
% correlation_plot_cell{2} = days to crit
% correlation_plot_cell{3} = amount learning
% correlation_plot_cell{4} = probe-over-probe similarity



datamtx_adj_problem = [rscore_predict_PreprobeProblem; rscore_unpredict_PreprobeProblem];
datamtx_last_probes = [rscore_predict_PreprobeFinalprobe; rscore_unpredict_PreprobeFinalprobe];
datamtx_last_problem = [rscore_predict_ProblemFinalprobe; rscore_unpredict_ProblemFinalprobe];
num_subjects = length(predict_subjs) + length(unpredict_subjs); 
num_sessions = length(preprobe_session_numbers);
subj_num_mtx = repmat((1:num_subjects)', 1, num_sessions); % subject matrix
session_num_mtx = repmat(1:num_sessions, num_subjects, 1); % stage matrix
group_mtx = [zeros(length(predict_subjs), num_sessions); ones(length(unpredict_subjs), num_sessions)];
    
%
% First day discrimination vs Overlap (preprobe and problem (mean of first day and last day))  
%model_str = 'PopOverlap~Discrim+(1|Subject)+(1|Session)';
%tbl = table(datamtx_adj_problem(:), first_day_discrimination(:), subj_num_mtx(:), session_num_mtx(:), 'VariableNames', {'PopOverlap', 'Discrim', 'Subject', 'Session'});

    model_str = 'PopOverlap~Discrim+(1|Subject)';
    tbl = table(datamtx_adj_problem(:), first_day_discrimination(:), subj_num_mtx(:), 'VariableNames', {'PopOverlap', 'Discrim', 'Subject'});


    tbl.Subject = categorical(tbl.Subject);
    lme_PreprobeProblem_corr = fitlme(tbl, model_str)
   
    figure(f1); subplot(1,4,4); hold on
    plot2var(reshape(first_day_discrimination(group_mtx==0), length(predict_subjs), 6), rscore_predict_PreprobeProblem, ...
                reshape(first_day_discrimination(group_mtx==1), length(unpredict_subjs), 6), rscore_unpredict_PreprobeProblem, 'ProblemPreprobe')
            xlim([-5 5]); ylim([-0.1 1.1]); yticks(0:.2:1)
            xlabel('First day discrimination (d`)')
            ylabel('r-score')
            hold on; plot(xlim, [0 0], 'k--'); plot([0 0], ylim, 'k--')
            var_name = 'rscore_PreprobeProblem'; 
            print([figure_path var_name], '-dpdf', '-painters', '-bestfit')

% training reactivation vs final reactivation(problem and final probe)
%model_str = 'PopOverlap~Discrim+(1|Subject)+(1|Session)';
%tbl = table(datamtx_last_problem(:), first_day_discrimination(:), subj_num_mtx(:), session_num_mtx(:), 'VariableNames', {'PopOverlap', 'Discrim', 'Subject', 'Session'});

    model_str = 'PopOverlap~Discrim+(1|Subject)';
    tbl = table(datamtx_last_problem(:), datamtx_adj_problem(:), subj_num_mtx(:), 'VariableNames', {'PopOverlap', 'Discrim', 'Subject'});
    
    tbl.Subject = categorical(tbl.Subject);
    lme_ProblemFinalprobe_corr = fitlme(tbl, model_str)

    figure(f3); subplot(1,4,4); hold on
    plot2var(rscore_predict_PreprobeProblem, rscore_predict_ProblemFinalprobe, ...
                rscore_unpredict_PreprobeProblem, rscore_unpredict_ProblemFinalprobe, 'ProblemFinalprobe')
            axis auto; xlim([-.1 1.1]); ylim([-0.1 1.1]) ; yticks(0:.2:1)
            xlabel('r-score during training')
            ylabel('r-score with final probe')
            hold on; plot(xlim, [0 0], 'k--'); plot([0 0], ylim, 'k--')
            var_name = 'rscore_ProblemFinalprobe'; 
            print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
            
            
% training reactivation vs final reactivation (preprobe and final probe)
%model_str = 'PopOverlap~Discrim+(1|Subject)+(1|Session)';
%tbl = table(datamtx_last_probes(:), first_day_discrimination(:), subj_num_mtx(:), session_num_mtx(:), 'VariableNames', {'PopOverlap', 'Discrim', 'Subject', 'Session'});

    model_str = 'FinalReactivation~TrainingReactivation+(1|Subject)';
    tbl = table(datamtx_last_probes(:), datamtx_adj_problem(:), subj_num_mtx(:), 'VariableNames', {'FinalReactivation', 'TrainingReactivation', 'Subject'});
    
    tbl.Subject = categorical(tbl.Subject);
    lme_PreprobeFinalprobe_corr = fitlme(tbl, model_str)

    figure(f2); subplot(1,4,4); hold on
    plot2var(rscore_predict_PreprobeProblem, rscore_predict_PreprobeFinalprobe, ...
                rscore_unpredict_PreprobeProblem, rscore_unpredict_PreprobeFinalprobe, 'PreprobeFinalprobe')
            axis auto; xlim([-.1 1.1]); ylim([-0.1 1.1]); yticks(0:.2:1) 
            xlabel('r-score during training')
            ylabel('r-score with final probe')
            hold on; plot(xlim, [0 0], 'k--'); plot([0 0], ylim, 'k--')
            var_name = 'rscore_PreprobeFinalprobe'; 
            print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
            


%}            
            
            
%% Excess overlap with final probe in terms of excess overlap during training
%{
% First day discrimination vs Overlap (preprobe and final probe)
model_str = 'PopOverlap~Discrim+(1|Subject)';

    tbl = table(datamtx_last_probes(:), datamtx_adj_problem(:), subj_num_mtx(:), 'VariableNames', {'PopOverlap', 'Discrim', 'Subject'});
    tbl.Subject = categorical(tbl.Subject);
    lme_PreprobeFinalprobe_vs_trainingOverlap = fitlme(tbl, model_str)

    plot2var(prop_predict_PreprobeProblem, prop_predict_PreprobeFinalprobe, ...
                prop_unpredict_PreprobeProblem, prop_unpredict_PreprobeFinalprobe, 'Preprobe-LastProbe overlap vs Preprobe-Problem overlap')
            axis([-.5 .5 -.3 .4]); 
            xlabel('First day discrimination (d`)')
            ylabel('Population overlap')
            hold on; plot(xlim, [0 0], 'k--'); plot([0 0], ylim, 'k--')
            var_name = 'CorrDot_OverlapPreprobeLastprobe_vs_trainingoverlap'; 
            print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
            

% First day discrimination vs Overlap (problem and final probe)
model_str = 'PopOverlap~Discrim+(1|Subject)';

    tbl = table(datamtx_last_problem(:), datamtx_adj_problem(:), subj_num_mtx(:), 'VariableNames', {'PopOverlap', 'Discrim', 'Subject'});
    tbl.Subject = categorical(tbl.Subject);
    lme_ProblemFinalprobe_vs_trainingOverlap = fitlme(tbl, model_str)

    plot2var(prop_predict_PreprobeProblem, prop_predict_ProblemFinalprobe, ...
                prop_unpredict_PreprobeProblem, prop_unpredict_ProblemFinalprobe, 'Problem-LastProbe overlap vs Preprobe-Problem overlap')
            axis([-.5 .5 -.3 .4]); 
            xlabel('First day discrimination (d`)')
            ylabel('Population overlap')
            hold on; plot(xlim, [0 0], 'k--'); plot([0 0], ylim, 'k--')
            var_name = 'CorrDot_OverlapProblemLastprobe_vs_trainingoverlap'; 
            print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
%}            
            
            


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


function out = colon_op(in)

out = in(:);

end

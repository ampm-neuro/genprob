green_color = [21 154 72]./255; blue_color =  [43 73 159]./255;
light_green_color = [162 201 129]./255; light_blue_color = [102 164 213]./255;
figure_path = 'C:\Users\ampm1\Documents\manuscripts\gentrain\manuscript_figures\fig04\';
predict_subjs = [{'651049m1'} {'658648m2'} {'691359m1'} {'690330m1'} {'691359m2'} {'696944m3'} {'699437m4'} {'699438m3'}];
unpredict_subjs = [{'687034m4'} {'683472m2'} {'699437m2'} {'699438m1'} {'683472m3'}];


%% cell overlay example
%
load('C:\Users\ampm1\Documents\MATLAB\generalization\cellreg_data\696944m3\cell_reg_696944m3.mat');
    % great sessions: 1&2, 4&5, 10&11, 16&17
for isesh = 1 % session 1 (and 2) was used in figure
    figure; hold on
    trace_footprints_sessions_overlap(cell_registered_struct.spatial_footprints_corrected([isesh isesh+1]), cell_registered_struct.cell_to_index_map(:,[isesh isesh+1]))
    legend('location', 'northeastoutside')
    legend({'session 2', 'session 1'})
    var_name = 'Overlapping_population_example_venn'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit'); close
    var_name = 'Overlapping_population_example_footprints'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
end
%}


%% cell origin color plot (waterfall)

% rainbow colors
num_colors = 20;
colors = hsv(num_colors);

% sessions chronological order index
session_chron_reorder_idx = [1 9 15 2 10 16 3 11 17 4 12 18 5 13 19 6 14 20 7 8];

% combine cell registration matrices (sessions 1:20)
%{
[super_crm_out_predict, subject_cell_crm_predict] = super_crm(predict_subjs);
super_crm_out_predict = treat_crm(super_crm_out_predict, session_chron_reorder_idx);
for isubj = 1:length(subject_cell_crm_predict)
    subject_cell_crm_predict{isubj} = treat_crm(subject_cell_crm_predict{isubj}, session_chron_reorder_idx);
end
[super_crm_out_unpredict, subject_cell_crm_unpredict] = super_crm(unpredict_subjs);
super_crm_out_unpredict = treat_crm(super_crm_out_unpredict, session_chron_reorder_idx);
for isubj = 1:length(subject_cell_crm_unpredict)
    subject_cell_crm_unpredict{isubj} = treat_crm(subject_cell_crm_unpredict{isubj}, session_chron_reorder_idx);
end
save('cell_regist_fig04.mat', 'super_crm_out_predict', 'super_crm_out_unpredict', 'subject_cell_crm_predict', 'subject_cell_crm_unpredict')
%}
load('cell_regist_fig04.mat', 'super_crm_out_predict', 'super_crm_out_unpredict')

% Rainbow waterfall predict
figure; subplot(1,3,1:2)
waterfall_plot(super_crm_out_predict, colors); 
title('Predict')

    % bar chart of counts
    subplot(1,3,3); 
    rainbow_bar(super_crm_out_predict, colors)
    ylim([0 500]); yticks([0 250 500])
    
    % save
    set(gcf, 'Position', [607   559   633   420])
    var_name = 'Rainbow_predict'; 
    print([figure_path var_name], '-dpdf', '-painters', '-bestfit')

% Rainbow waterfall unpredict
    figure; subplot(1,3,1:2)
    waterfall_plot(super_crm_out_unpredict, colors); 
    title('Unpredict')

    % bar chart of counts
    subplot(1,3,3); 
    rainbow_bar(super_crm_out_unpredict, colors)
    ylim([0 500]); yticks([0 250 500])
    
    % save
    set(gcf, 'Position', [607   559   633   420])
    var_name = 'Rainbow_unpredict'; 
    print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
%}    

%% proportion of overlap matrices

load('cell_regist_fig04.mat', 'subject_cell_crm_predict', 'subject_cell_crm_unpredict')
session_numbers = [1:1:19];

% preallcate sessions x sessions x subjects
prop_predict = nan(length(session_numbers), length(session_numbers), length(predict_subjs));
prop_unpredict = nan(length(session_numbers), length(session_numbers), length(unpredict_subjs));

% iterate through sessions
for isesh1 = 1:length(session_numbers)
    for isesh2 = 1:length(session_numbers)
        
        % iterate through predict subjects
        for isubj_predict = 1:length(predict_subjs)
            % load crm
            crm = subject_cell_crm_predict{isubj_predict};
            % number unique cells
            unq_ct = sum(~isnan(crm(:,session_numbers(isesh1))) | ~isnan(crm(:,session_numbers(isesh2))));
            % active in both
            common_ct = sum(~isnan(crm(:,session_numbers(isesh1))) & ~isnan(crm(:,session_numbers(isesh2))));
            % load 
            prop_predict(isesh1, isesh2, isubj_predict) = common_ct/unq_ct;
        end
        
        % iterate through unpredict subjects
        for isubj_unpredict = 1:length(unpredict_subjs)
            % load crm
            crm = subject_cell_crm_unpredict{isubj_unpredict};
            % number unique cells
            unq_ct = sum(~isnan(crm(:,session_numbers(isesh1))) | ~isnan(crm(:,session_numbers(isesh2))));
            % active in both
            common_ct = sum(~isnan(crm(:,session_numbers(isesh1))) & ~isnan(crm(:,session_numbers(isesh2))));
            % load 
            prop_unpredict(isesh1, isesh2, isubj_unpredict) = common_ct/unq_ct;
        end
        
    end
end
    
figure;
subplot(1,2,1); imagesc(nanmean(prop_predict,3))
axis square; caxis([0 .4]); set(gca,'TickLength',[0, 0]); colorbar
xticks(1:19); yticks(1:19); title('Predictable training')
subplot(1,2,2); imagesc(nanmean(prop_unpredict,3))
axis square; caxis([0 .4]); set(gca,'TickLength',[0, 0]); colorbar
xticks(1:19); yticks(1:19); title('Unpredictable training')
var_name = 'Mean_overlap_matrices'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')


save('overlap_raw', 'prop_predict', 'prop_unpredict')




%% correct for time

% time between imaging session (in days)
%{
days_mtx_predict = [];
for isubj_predict = 1:length(predict_subjs)
    days_mtx_predict = cat(3, days_mtx_predict, session_days(predict_subjs{isubj_predict}));
end
days_mtx_unpredict = [];
for isubj_unpredict = 1:length(unpredict_subjs)
    days_mtx_unpredict = cat(3, days_mtx_unpredict, session_days(unpredict_subjs{isubj_unpredict}));
end
figure; subplot(1,2,1); imagesc(nanmean(days_mtx_predict,3)); caxis([0 60]);...
 subplot(1,2,2); imagesc(nanmean(days_mtx_unpredict,3)); caxis([0 60]);
save('days_between.mat', 'days_mtx_predict', 'days_mtx_unpredict')
%}
load('overlap_raw', 'prop_predict', 'prop_unpredict')
load('days_between.mat', 'days_mtx_predict', 'days_mtx_unpredict')

% remove identity diaganol
prop_predict(prop_predict>0.99) = nan; 
prop_unpredict(prop_unpredict>0.99) = nan; 
days_mtx_predict = days_mtx_predict(1:19,1:19,:);
    days_mtx_predict(days_mtx_predict==0) = nan;
days_mtx_unpredict = days_mtx_unpredict(1:19,1:19,:);
    days_mtx_unpredict(days_mtx_unpredict==0) = nan;

max_x = 100;
% combine
%all_overlap = [prop_predict(:); prop_unpredict(:)]; 
    all_overlap = [prop_unpredict(:)]; 
%all_days = [days_mtx_predict(:);days_mtx_unpredict(:)];
    all_days = [days_mtx_unpredict(:)];
all_overlap = all_overlap(~isnan(all_days) & all_days<=max_x);
all_days = all_days(~isnan(all_days) & all_days<=max_x);
    
% model time
ibins = [0:0.25:5 5.5:.5:30 31:1:40 42:2:60 65:5:100];
ibin_cell = cell(1,length(ibins)-1);
for ibin = 1:length(ibin_cell)
    ibin_cell{ibin} = all_overlap(all_days>ibins(ibin) & all_days<=ibins(ibin+1));
end
[~, coefEsts, modelFun] = ampm_normal_fit(all_days(~isnan(all_days) & all_days<max_x), all_overlap(~isnan(all_days) & all_days<max_x), [0 1 0 150]);

% correct population overlap matrices
prop_predict = prop_predict - modelFun(coefEsts, days_mtx_predict);
prop_unpredict = prop_unpredict - modelFun(coefEsts, days_mtx_unpredict);



% model comparing overlap matrices
datamtx = cat(3, prop_predict, prop_unpredict);
for isesh1 = 1:size(prop_predict,1)
    datamtx(isesh1, isesh1, :) = nan;
end

num_subjects = length(predict_subjs)+length(unpredict_subjs); 
num_sessions1 = 19;
num_sessions2 = 19;
subj_num_mtx = repmat(permute(1:num_subjects, [1 3 2]), [num_sessions1,num_sessions2,1]);
session_num_mtx1 = repmat((1:num_sessions1)', [1,num_sessions2,num_subjects]);
session_num_mtx2 = repmat(1:num_sessions1, [num_sessions1,1,num_subjects]);
group_mtx = cat(3,zeros(num_sessions1,num_sessions2,length(predict_subjs)), ones(num_sessions1,num_sessions2,length(unpredict_subjs)));

nnan_idx = ~isnan(datamtx);
datamtx = datamtx(nnan_idx);
subj_num_mtx = subj_num_mtx(nnan_idx);
session_num_mtx1 = session_num_mtx1(nnan_idx);
session_num_mtx2 = session_num_mtx2(nnan_idx);
group_mtx = group_mtx(nnan_idx);

model_str = 'MemorySimilarity~Session1*Session2*Group+(1|Subject)';
tbl = table(datamtx, session_num_mtx1, session_num_mtx2, group_mtx, subj_num_mtx, 'VariableNames', {'MemorySimilarity', 'Session1', 'Session2', 'Group', 'Subject'});
tbl.Subject = categorical(tbl.Subject);
lme_PreprobeFinalprobe = fitlme(tbl, model_str)

    % follow up for predict group
    model_str = 'MemorySimilarity~Session1*Session2+(1|Subject)';
    tbl = table(datamtx(group_mtx==0), session_num_mtx1(group_mtx==0), session_num_mtx2(group_mtx==0), subj_num_mtx(group_mtx==0), 'VariableNames', {'MemorySimilarity', 'Session1', 'Session2', 'Subject'});
    tbl.Subject = categorical(tbl.Subject);
    lme_PreprobeFinalprobe_predict = fitlme(tbl, model_str)
    
    % follow up for unpredict group
    model_str = 'MemorySimilarity~Session1*Session2+(1|Subject)';
    tbl = table(datamtx(group_mtx==1), session_num_mtx1(group_mtx==1), session_num_mtx2(group_mtx==1), subj_num_mtx(group_mtx==1), 'VariableNames', {'MemorySimilarity', 'Session1', 'Session2', 'Subject'});
    tbl.Subject = categorical(tbl.Subject);
    lme_PreprobeFinalprobe_unpredict = fitlme(tbl, model_str)









%% Overlap comparisons (errorbars)
%{
preprobe_session_numbers = 1:3:16;
firstday_session_numbers = 2:3:17;
lastday_session_numbers = 3:3:18;
postprobe_session_numbers = 4:3:19;

% preallcate subjects x probes
prop_predict_PreprobeProblem = nan(length(predict_subjs), length(preprobe_session_numbers));
prop_predict_ProblemFinalprobe = nan(length(predict_subjs), length(preprobe_session_numbers));
prop_unpredict_PreprobeProblem = nan(length(unpredict_subjs), length(preprobe_session_numbers));
prop_unpredict_ProblemFinalprobe = nan(length(unpredict_subjs), length(preprobe_session_numbers));
prop_predict_PreprobeFinalprobe = nan(length(predict_subjs), length(preprobe_session_numbers));
prop_unpredict_PreprobeFinalprobe = nan(length(unpredict_subjs), length(preprobe_session_numbers));

% iterate through sessions
for isesh = 1:length(preprobe_session_numbers)
    
    % adjacent probe sessions
    %
        % iterate through predict subjects
        for isubj_predict = 1:length(predict_subjs)
            prop_predict_PreprobeProblem(isubj_predict, isesh) = ...
                mean(prop_predict(preprobe_session_numbers(isesh), [firstday_session_numbers(isesh) lastday_session_numbers(isesh)], isubj_predict));
        end
        % iterate through unpredict subjects
        for isubj_unpredict = 1:length(unpredict_subjs)
            prop_unpredict_PreprobeProblem(isubj_unpredict, isesh) = ...
                mean([prop_unpredict(preprobe_session_numbers(isesh), firstday_session_numbers(isesh), isubj_unpredict);...
                prop_unpredict(preprobe_session_numbers(isesh), lastday_session_numbers(isesh), isubj_unpredict);]);
        end
        
    % compare all problems with last probe
    %
        % iterate through predict subjects
        for isubj_predict = 1:length(predict_subjs)
            prop_predict_ProblemFinalprobe(isubj_predict, isesh) = ...
                mean([prop_predict(firstday_session_numbers(isesh), end, isubj_predict);... 
                      prop_predict(lastday_session_numbers(isesh), end, isubj_predict)]); % compared to probe 6
        end
        % iterate through unpredict subjects
        for isubj_unpredict = 1:length(unpredict_subjs)
            prop_unpredict_ProblemFinalprobe(isubj_unpredict, isesh) = ...
                mean([prop_unpredict(firstday_session_numbers(isesh), end, isubj_unpredict);...
                      prop_unpredict(lastday_session_numbers(isesh), end, isubj_unpredict)]);% compared to probe 6
        end
        
    % compare all probes with last probe
    %
        % iterate through predict subjects
        for isubj_predict = 1:length(predict_subjs)
            prop_predict_PreprobeFinalprobe(isubj_predict, isesh) = prop_predict(preprobe_session_numbers(isesh), postprobe_session_numbers(end), isubj_predict);
        end
        % iterate through unpredict subjects
        for isubj_unpredict = 1:length(unpredict_subjs)
            prop_unpredict_PreprobeFinalprobe(isubj_unpredict, isesh) = prop_unpredict(preprobe_session_numbers(isesh), postprobe_session_numbers(end), isubj_unpredict);
        end
end

% combine data for later mixed models
datamtx_adj_PreprobeProblem = [prop_predict_PreprobeProblem; prop_unpredict_PreprobeProblem];
datamtx_last_PreprobeFinalprobe = [prop_predict_PreprobeFinalprobe; prop_unpredict_PreprobeFinalprobe];
datamtx_last_ProblemFinalprobe = [prop_predict_ProblemFinalprobe; prop_unpredict_ProblemFinalprobe];

% plot preprobe vs problem
%
figure; hold on; 
errorbar_plot(num2cell(prop_predict_PreprobeProblem,1), 1, [], light_green_color, green_color)
errorbar_plot(num2cell(prop_unpredict_PreprobeProblem,1), 1, [], light_blue_color, blue_color)
xticks(1:size(prop_predict_PreprobeProblem,2)); xlabel('Probe number'); xticklabels(0:size(prop_predict_PreprobeFinalprobe,2)-1);
ylim([-.3 .4]); yticks([-1:.1:1]); ylabel('Population overlap')
axis square; title('preprobe and problem')
hold on; plot(xlim, [0 0], 'k--');

    % mixed model
    num_subjects = length(predict_subjs) + length(unpredict_subjs); 
    num_sessions = length(preprobe_session_numbers);
    subj_num_mtx = repmat((1:num_subjects)', 1, num_sessions); % subject matrix
    session_num_mtx = repmat(1:num_sessions, num_subjects, 1); % stage matrix
    group_mtx = [zeros(length(predict_subjs), num_sessions); ones(length(unpredict_subjs), num_sessions)];
    model_str = 'PopOverlap~Session*Group+(1|Subject)';
    tbl = table(datamtx_adj_PreprobeProblem(:), session_num_mtx(:), group_mtx(:), subj_num_mtx(:), 'VariableNames', {'PopOverlap', 'Session', 'Group', 'Subject'});
    tbl.Subject = categorical(tbl.Subject);
    lme_PreprobeProblem = fitlme(tbl, model_str);

        % one way anova
        rmtbl_adj_PreprobeProblem_predict = simple_mixed_anova(prop_predict_PreprobeProblem);
        rmtbl_adj_PreprobeProblem_unpredict = simple_mixed_anova(prop_unpredict_PreprobeProblem);
        
    % ttest compare to zero
    pval_predict_PreprobeProblem_0 = nan(size(prop_predict_PreprobeProblem,2),1); 
        stats_predict_PreprobeProblem_0 = cell(size(prop_predict_PreprobeProblem,2),1);
    pval_unpredict_PreprobeProblem_0 = nan(size(prop_unpredict_PreprobeProblem,2),1); 
        stats_unpredict_PreprobeProblem_0 = cell(size(prop_unpredict_PreprobeProblem,2),1);
    for icomp = 1:size(prop_predict_PreprobeProblem,2)
        [~,pval_predict_PreprobeProblem_0(icomp),~,stats_predict_PreprobeProblem_0{icomp}] = ttest(prop_predict_PreprobeProblem(:,icomp));
        [~,pval_unpredict_PreprobeProblem_0(icomp),~,stats_unpredict_PreprobeProblem_0{icomp}] = ttest(prop_unpredict_PreprobeProblem(:,icomp));
    end

    % ttest stats within groups
    pval_adj_predict_firstlast = nan(size(prop_predict_PreprobeProblem,2),1); stats_adj_predict_firstlast = cell(size(prop_predict_PreprobeProblem,2),1);
    pval_adj_unpredict_firstlast = nan(size(prop_unpredict_PreprobeProblem,2),1); stats_adj_unpredict_firstlast = cell(size(prop_unpredict_PreprobeProblem,2),1);
    for icomp = 2:size(prop_predict_PreprobeProblem,2)
        [~,pval_adj_predict_firstlast(icomp),~,stats_adj_predict_firstlast{icomp}] = ttest(prop_predict_PreprobeProblem(:,1), prop_predict_PreprobeProblem(:,icomp));
        [~,pval_adj_unpredict_firstlast(icomp),~,stats_adj_unpredict_firstlast{icomp}] = ttest(prop_unpredict_PreprobeProblem(:,1), prop_unpredict_PreprobeProblem(:,icomp));
    end

    % ttest stats between groups
    pval_adj_comp_firstlast = nan(size(prop_unpredict_PreprobeProblem,2),1); stats_adj_comppval_adj_comp_firstlast = cell(size(prop_unpredict_PreprobeProblem,2),1);
    for icomp = 1:size(prop_predict_PreprobeProblem,2)
        [~,pval_adj_comp_firstlast(icomp),~,stats_adj_comppval_adj_comp_firstlast{icomp}] = ttest2(prop_predict_PreprobeProblem(:,icomp), prop_unpredict_PreprobeProblem(:,icomp));
        pval_text(pval_adj_comp_firstlast(icomp), icomp-.4, 0.35)
    end
    
    var_name = 'overlap_PreprobeProblem'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')      

% plot preprobe vs probe 6
%
figure; hold on; 
errorbar_plot(num2cell(prop_predict_PreprobeFinalprobe,1), 1, [], light_green_color, green_color)
errorbar_plot(num2cell(prop_unpredict_PreprobeFinalprobe,1), 1, [], light_blue_color, blue_color)
xticks(1:size(prop_predict_PreprobeFinalprobe,2)); xlabel('Probe number'); xticklabels(0:size(prop_predict_PreprobeFinalprobe,2)-1);
ylim([-.3 .4]); yticks([-1:.1:1]); ylabel('Population overlap')
axis square; title('probes to last probe')
hold on; plot(xlim, [0 0], 'k--');

    % mixed model 
    num_subjects = length(predict_subjs) + length(unpredict_subjs); 
    num_sessions = length(preprobe_session_numbers);
    subj_num_mtx = repmat((1:num_subjects)', 1, num_sessions); % subject matrix
    session_num_mtx = repmat(1:num_sessions, num_subjects, 1); % stage matrix
    group_mtx = [zeros(length(predict_subjs), num_sessions); ones(length(unpredict_subjs), num_sessions)];
    model_str = 'PopOverlap~Session*Group+(1|Subject)';
    tbl = table(datamtx_last_PreprobeFinalprobe(:), session_num_mtx(:), group_mtx(:), subj_num_mtx(:), 'VariableNames', {'PopOverlap', 'Session', 'Group', 'Subject'});
    tbl.Subject = categorical(tbl.Subject);
    lme_PreprobeFinalprobe = fitlme(tbl, model_str)
    
        % one way anova
        rmtbl_adj_PreprobeFinalprobe_predict = simple_mixed_anova(prop_predict_PreprobeFinalprobe);
        rmtbl_adj_PreprobeFinalprobe_unpredict = simple_mixed_anova(prop_unpredict_PreprobeFinalprobe);

    % ttest compare to zero
    pval_last_predict_0 = nan(size(prop_predict_PreprobeFinalprobe,2),1); stats_last_predict_0 = cell(size(prop_predict_PreprobeFinalprobe,2),1);
    pval_last_unpredict_0 = nan(size(prop_unpredict_PreprobeFinalprobe,2),1); stats_last_unpredict_0 = cell(size(prop_unpredict_PreprobeFinalprobe,2),1);
    for icomp = 1:size(prop_predict_PreprobeFinalprobe,2)
        [~,pval_last_predict_0(icomp),~,stats_last_predict_0{icomp}] = ttest(prop_predict_PreprobeFinalprobe(:,icomp));
        [~,pval_last_unpredict_0(icomp),~,stats_last_unpredict_0{icomp}] = ttest(prop_unpredict_PreprobeFinalprobe(:,icomp));
    end

    % ttest stats within groups
    pval_last_predict = nan(size(prop_predict_PreprobeFinalprobe,2),1); stats_last_predict = cell(size(prop_predict_PreprobeFinalprobe,2),1);
    pval_last_unpredict = nan(size(prop_unpredict_PreprobeFinalprobe,2),1); stats_last_unpredict = cell(size(prop_unpredict_PreprobeFinalprobe,2),1);
    for icomp = 2:size(prop_predict_PreprobeFinalprobe,2)
        [~,pval_last_predict(icomp),~,stats_last_predict{icomp}] = ttest(prop_predict_PreprobeFinalprobe(:,1), prop_predict_PreprobeFinalprobe(:,icomp));
        [~,pval_last_unpredict(icomp),~,stats_last_unpredict{icomp}] = ttest(prop_unpredict_PreprobeFinalprobe(:,1), prop_unpredict_PreprobeFinalprobe(:,icomp));
    end

    % ttest stats between group
    pval_last_comp = nan(size(prop_unpredict_PreprobeFinalprobe,2),1); stats_last_comp = cell(size(prop_unpredict_PreprobeFinalprobe,2),1);
    for icomp = 1:size(prop_predict_PreprobeFinalprobe,2)
        [~,pval_last_comp(icomp),~,stats_last_comp{icomp}] = ttest2(prop_predict_PreprobeFinalprobe(:,icomp), prop_unpredict_PreprobeFinalprobe(:,icomp));
        pval_text(pval_last_comp(icomp), icomp-.4, 0.35)
    end
    
    var_name = 'overlap_PreprobeFinalprobe'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
    
    
% plot problem vs probe 6
%
figure; hold on; 
errorbar_plot(num2cell(prop_predict_ProblemFinalprobe,1), 1, [], light_green_color, green_color)
errorbar_plot(num2cell(prop_unpredict_ProblemFinalprobe,1), 1, [], light_blue_color, blue_color)
xticks(1:size(prop_predict_ProblemFinalprobe,2)); xlabel('Problem number'); xticklabels(0:size(prop_predict_ProblemFinalprobe,2)-1);
ylim([-.3 .4]); yticks([-1:.1:1]); ylabel('Population overlap')
axis square; title('problem and last probe')
hold on; plot(xlim, [0 0], 'k--');
 
    % mixed model
    num_subjects = length(predict_subjs) + length(unpredict_subjs); 
    num_sessions = length(preprobe_session_numbers);
    subj_num_mtx = repmat((1:num_subjects)', 1, num_sessions); % subject matrix
    session_num_mtx = repmat(1:num_sessions, num_subjects, 1); % stage matrix
    group_mtx = [zeros(length(predict_subjs), num_sessions); ones(length(unpredict_subjs), num_sessions)];
    model_str = 'PopOverlap~Session*Group+(1|Subject)';
    tbl = table(datamtx_last_ProblemFinalprobe(:), session_num_mtx(:), group_mtx(:), subj_num_mtx(:), 'VariableNames', {'PopOverlap', 'Session', 'Group', 'Subject'});
    tbl.Subject = categorical(tbl.Subject);
    lme_ProblemFinalprobe = fitlme(tbl, model_str)

        % one way anova
        rmtbl_adj_ProblemFinalprobe_predict = simple_mixed_anova(prop_predict_ProblemFinalprobe);
        rmtbl_adj_ProblemFinalprobe_unpredict = simple_mixed_anova(prop_unpredict_ProblemFinalprobe);

    % ttest compare to zero
    pval_last_predict_firstlast_0 = nan(size(prop_predict_ProblemFinalprobe,2),1); 
        stats_last_predict_firstlast_0 = cell(size(prop_predict_ProblemFinalprobe,2),1);
    pval_last_unpredict_firstlast_0 = nan(size(prop_unpredict_ProblemFinalprobe,2),1); 
        stats_unpredict_firstlast_0 = cell(size(prop_unpredict_ProblemFinalprobe,2),1);
    for icomp = 1:size(prop_predict_ProblemFinalprobe,2)
        [~,pval_last_predict_firstlast_0(icomp),~,stats_last_predict_firstlast_0{icomp}] = ttest(prop_predict_ProblemFinalprobe(:,icomp));
        [~,pval_last_unpredict_firstlast_0(icomp),~,stats_unpredict_firstlast_0{icomp}] = ttest(prop_unpredict_ProblemFinalprobe(:,icomp));
    end

    % ttest stats within groups
    pval_last_predict_firstlast = nan(size(prop_predict_ProblemFinalprobe,2),1); 
        stats_last_predict_firstlast = cell(size(prop_predict_ProblemFinalprobe,2),1);
    pval_last_unpredict_firstlast = nan(size(prop_unpredict_ProblemFinalprobe,2),1); 
        stats_last_unpredict_firstlast = cell(size(prop_unpredict_ProblemFinalprobe,2),1);
    for icomp = 2:size(prop_predict_ProblemFinalprobe,2)
        [~,pval_last_predict_firstlast(icomp),~,stats_last_predict_firstlast{icomp}] = ttest(prop_predict_ProblemFinalprobe(:,1), prop_predict_ProblemFinalprobe(:,icomp));
        [~,pval_last_unpredict_firstlast(icomp),~,stats_last_unpredict_firstlast{icomp}] = ttest(prop_unpredict_ProblemFinalprobe(:,1), prop_unpredict_ProblemFinalprobe(:,icomp));
    end

    % ttest stats between group
    pval_last_comp_firstlast = nan(size(prop_unpredict_ProblemFinalprobe,2),1); stats_last_comp_firstlast = cell(size(prop_unpredict_ProblemFinalprobe,2),1);
    for icomp = 1:size(prop_predict_ProblemFinalprobe,2)
        [~,pval_last_comp_firstlast(icomp),~,stats_last_comp_firstlast{icomp}] = ttest2(prop_predict_ProblemFinalprobe(:,icomp), prop_unpredict_ProblemFinalprobe(:,icomp));
        pval_text(pval_last_comp_firstlast(icomp), icomp-.4, 0.35)
    end
    
    var_name = 'overlap_ProblemFinalprobe'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')



%% Overlap comparisons (correlations)
%
%[~, ~, diff_waits_problem_preprobe_predict, diff_waits_problem_postprobe_predict, correlation_plot_cell_predict] = ALL_probe_vs_behavior_fit('train_mevar_imaging_hpc');
%for ifig = 1:34; close; end
%[~, ~, diff_waits_problem_preprobe_unpredict, diff_waits_problem_postprobe_unpredict, correlation_plot_cell_unpredict] = ALL_probe_vs_behavior_fit('train_hivar_imaging_hpc');
%for ifig = 1:34; close; end
load('sfigure02_2', 'diff_waits_problem_preprobe_predict', 'diff_waits_problem_preprobe_unpredict', 'correlation_plot_cell_predict', 'correlation_plot_cell_unpredict')

accuracy_of_prediction = [diff_waits_problem_preprobe_predict; diff_waits_problem_preprobe_unpredict];
first_day_discrimination = [correlation_plot_cell_predict{1}; correlation_plot_cell_unpredict{1}];
datamtx_adj_problem = [prop_predict_PreprobeProblem; prop_unpredict_PreprobeProblem];
datamtx_last_problem = [prop_predict_ProblemFinalprobe; prop_unpredict_ProblemFinalprobe];
datamtx_last_probes = [prop_predict_PreprobeFinalprobe; prop_unpredict_PreprobeFinalprobe];


% First day discrimination vs Overlap (preprobe and problem (mean of first day and last day))  
model_str = 'PopOverlap~Discrim+(1|Subject)';

    tbl = table(datamtx_adj_problem(:), first_day_discrimination(:), subj_num_mtx(:), 'VariableNames', {'PopOverlap', 'Discrim', 'Subject'});
    tbl.Subject = categorical(tbl.Subject);
    lme_PreprobeProblem = fitlme(tbl, model_str)
   
    plot2var(correlation_plot_cell_predict{1}, prop_predict_PreprobeProblem, ...
                correlation_plot_cell_unpredict{1}, prop_unpredict_PreprobeProblem, 'First day discrimination and population overlap between preprobe and problem sessions')
            axis([-5 5 -.3 .4]); 
            xlabel('First day discrimination (d`)')
            ylabel('Population overlap')
            hold on; plot(xlim, [0 0], 'k--'); plot([0 0], ylim, 'k--')
            var_name = 'FirstdayDiscrim_OverlapPreprobeFirstLastproblem'; 
            print([figure_path var_name], '-dpdf', '-painters', '-bestfit')

            disp('HERE')
            
% First day discrimination vs Overlap (preprobe and final probe)
model_str = 'PopOverlap~Discrim+(1|Subject)';

    tbl = table(datamtx_last_probes(:), first_day_discrimination(:), subj_num_mtx(:), 'VariableNames', {'PopOverlap', 'Discrim', 'Subject'});
    tbl.Subject = categorical(tbl.Subject);
    lme_PreprobeFinalprobe = fitlme(tbl, model_str)

    plot2var(prop_predict_PreprobeProblem, prop_predict_PreprobeFinalprobe, ...
                prop_predict_PreprobeProblem, prop_unpredict_PreprobeFinalprobe, 'First day discrimination and population overlap between preprobe and last probe')
            axis([-5 5 -.3 .4]); 
            xlabel('First day discrimination (d`)')
            ylabel('Population overlap')
            hold on; plot(xlim, [0 0], 'k--'); plot([0 0], ylim, 'k--')
            var_name = 'FirstdayDiscrim_OverlapPreprobeLastprobe'; 
            print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
            

% First day discrimination vs Overlap (problem and final probe)
model_str = 'PopOverlap~Discrim+(1|Subject)';

    tbl = table(datamtx_last_problem(:), first_day_discrimination(:), subj_num_mtx(:), 'VariableNames', {'PopOverlap', 'Discrim', 'Subject'});
    tbl.Subject = categorical(tbl.Subject);
    lme_ProblemFinalprobe = fitlme(tbl, model_str)

    plot2var(prop_predict_PreprobeProblem, prop_predict_ProblemFinalprobe, ...
                prop_predict_PreprobeProblem, prop_unpredict_ProblemFinalprobe, 'First day discrimination and population overlap between problem sessions and last probe')
            axis([-5 5 -.3 .4]); 
            xlabel('First day discrimination (d`)')
            ylabel('Population overlap')
            hold on; plot(xlim, [0 0], 'k--'); plot([0 0], ylim, 'k--')
            var_name = 'FirstdayDiscrim_OverlapFirstLastproblemLastprobe'; 
            print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
           
            
            
%% Excess overlap with final probe in terms of excess overlap during training
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

% sort combined crms into one psuedo subject matrix
function crm = treat_crm(crm, session_chron_reorder_idx)

    crm = crm(:, session_chron_reorder_idx);
    crm = crm(:,1:19);
    crm(crm==0) = nan;
    crm(sum(~isnan(crm),2)==0,:) = [];

    % label cells by origin session
    for icell = 1:size(crm,1)
        crm(icell, ~isnan(crm(icell,:))) = find(~isnan(crm(icell,:)),1,'first')-1;
    end
    
    % sort by first session
    crm_new = nan(size(crm));
    ct = 0;
    for isesh = 1:size(crm,2)
        crm_new(ct+1:ct+sum(crm(:,isesh)==isesh-1), :) = crm(crm(:,isesh)==isesh-1,:);
        ct = ct+sum(crm(:,isesh)==isesh-1);
    end
    crm = crm_new;

    % sort within a session by reactivation probability
    for isesh = 1:size(crm,2)

        % just cells that originated during this session
        crm_sesh = crm(crm(:,isesh)==isesh-1,:);

        % reorder by number of reactivations
        [~, sort_idx] = sort(sum(crm_sesh==isesh-1,2), 'descend');
        crm(crm(:,isesh)==isesh-1,:) = crm_sesh(sort_idx,:);
    end

end

% waterfall
function waterfall_plot(crm, colors)
    ampm_pcolor(crm')
    yticks((1:19)+.5); yticklabels(1:19)
    xticks([1 500:500:size(crm,1) size(crm,1)])
    xlim([.5 size(crm,1)+.5])
    ylabel('Session'); xlabel('Neuron')
    colormap(colors(1:19,:))
    box off
end

% bar chart
function rainbow_bar(crm, colors)
    unq_sesh_ids = unique(crm(~isnan(crm)));
    bar_input = nan(size(crm,2), length(unq_sesh_ids));
    for isesh = 1:size(crm,2)
        bar_input(isesh, :) = histcounts(crm(:, isesh), [unq_sesh_ids; length(unq_sesh_ids)]);
    end
    
    ba = bar(bar_input, 'stacked', 'FaceColor', 'flat');
    for iba_color = 1:size(crm,2)
        ba(iba_color).CData = colors(iba_color,:);
    end
    
    set(gca,'TickLength',[0, 0]); box off;
    xticks([])
    xlim([0.5 size(crm,2)+0.5])
    ylabel('Neuron count')
    view([-90 90])
    set(gca, 'YDir','reverse','XDir','reverse')
end

% get number of days between sessions
function days_mtx = session_days(subject_id)

    % all subject sessions
    subj_sessions = [get_file_paths_targeted('E:\two_tone\train_mevar_imaging_hpc', ['Subject ' subject_id]); ...
        get_file_paths_targeted('E:\two_tone\train_hivar_imaging_hpc', ['Subject ' subject_id])];
    
    % get all imaging sessions in order
    imaging_files = [];
    for igen = 7:12    
        
        % gen files
        if igen < 10
            problem_sessions = subj_sessions(contains(subj_sessions, ['gen0' num2str(igen)]));
        else
            problem_sessions = subj_sessions(contains(subj_sessions, ['gen' num2str(igen)]));
        end
        
        % first
        imaging_files = [imaging_files; problem_sessions(1)];
        
        % last
        imaging_files = [imaging_files; problem_sessions(end)];
            
    end
    
    % plus probes
    imaging_files = [imaging_files; subj_sessions(contains(subj_sessions, 'probe'))];

    % find date of every file
    file_dates = [];
    for ifile = 1:size(imaging_files,1)
        
        % string bounds
        
            % date
            date_start = strfind(imaging_files{ifile},'!')+1;
            date_end = strfind(imaging_files{ifile},'_')-1;
            date_end = min(date_end(date_end>date_start));
            file_date = imaging_files{ifile}(date_start : date_end);
        
            % hour
            hour_start = strfind(imaging_files{ifile},'_')+1;
            hour_start = min(hour_start(hour_start>date_start));
            hour_end = strfind(imaging_files{ifile},'h')-1;
            hour_end = min(hour_end(hour_end>date_start));
            file_hour = imaging_files{ifile}(hour_start : hour_end);
            
            % minute
            min_start = strfind(imaging_files{ifile},'h')+1;
            min_start = min(min_start(min_start>date_start));
            min_end = strfind(imaging_files{ifile},'m')-1;
            min_end = min(min_end(min_end>date_start));
            file_min = imaging_files{ifile}(min_start : min_end);
        
        % file date
        file_dates = [file_dates; datetime([file_date ' ' file_hour ':' file_min ])];
            
    end
    
    % sort
    [file_dates_sorted, chron_idx] = sort(file_dates);
    imaging_files = imaging_files(chron_idx);
    
    % compute hours to each session
    days_mtx = nan(size(file_dates_sorted,1));
    for idate1 = 1:size(file_dates_sorted,1)
        for idate2 = 1:size(file_dates_sorted,1)
            days_mtx(idate1, idate2) = abs(hours(file_dates_sorted(idate1)-file_dates_sorted(idate2))/24);
        end
    end

end

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
    axis square
end
%}
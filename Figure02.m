%% load data

green_color = [21 154 72]./255; blue_color =  [43 73 159]./255;
light_green_color = [162 201 129]./255; light_blue_color = [102 164 213]./255;
figure_path = 'C:\Users\ampm1\Documents\manuscripts\gentrain\manuscript_figures\fig02\';




%% all probe curves
figure;

% predict
probe_fps_predict = [get_file_paths_targeted('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar', 'probe')];            
for iprobe = 1:7
    
    if iprobe<7
        probe_paths = probe_fps_predict(contains(probe_fps_predict, ['preprobe_0' num2str(iprobe)]));
    else
        probe_paths = probe_fps_predict(contains(probe_fps_predict, 'postprobe_01'));
    end
    subplot(2,7,iprobe)
    plot_allprobes(probe_paths); 
    rich_bounds_prob('mevar', iprobe);
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
probe_fps_predict = [get_file_paths_targeted('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_hivar', 'probe')];            
for iprobe = 1:7
    
    if iprobe<7
        probe_paths = probe_fps_predict(contains(probe_fps_predict, ['preprobe_0' num2str(iprobe)]));
    else
        probe_paths = probe_fps_predict(contains(probe_fps_predict, 'postprobe_01'));
    end
    subplot(2,7,iprobe+7)
    plot_allprobes(probe_paths); 
    rich_bounds_prob('hivar', iprobe);
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
%var_name = 'probe_curves'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')



%% correlation matrices and associated line plots
[probe_waits_predict, ~, cm_17_predict, offdiag_cell_predict, comp2one_cell_predict, comp2seven_cell_predict] = probe_corr_mtx('train_mevar');
%var_name = 'ProbeBeh_corm_predict'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')

    % compare every matrix item to 0
    pval_mtx_predict = nan(size(cm_17_predict,1), size(cm_17_predict,2));
    for icomp1 = 1:size(cm_17_predict,1)
        for icomp2 = 1:size(cm_17_predict,2)
            [~,pval_mtx_predict(icomp1,icomp2)] = ttest(cm_17_predict(icomp1, icomp2, :));
        end
    end

[probe_waits_unpredict, ~, cm_17_unpredict, offdiag_cell_unpredict, comp2one_cell_unpredict, comp2seven_cell_unpredict] = probe_corr_mtx('train_hivar');
%var_name = 'ProbeBeh_corm_unpredict'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')

    % compare every matrix item to 0
    pval_mtx_unpredict = nan(size(cm_17_unpredict,1), size(cm_17_unpredict,2));
    for icomp1 = 1:size(cm_17_unpredict,1)
        for icomp2 = 1:size(cm_17_unpredict,2)
            [~,pval_mtx_unpredict(icomp1,icomp2)] = ttest(cm_17_unpredict(icomp1, icomp2, :));
        end
    end

% model comparing correlation matrices
cm_17_predict_noduplicates = cm_17_predict;
cm_17_unpredict_noduplicates = cm_17_unpredict;
for isesh1 = 1:7
    cm_17_predict_noduplicates(isesh1, isesh1, :) = nan;
    cm_17_unpredict_noduplicates(isesh1, isesh1, :) = nan;
end
datamtx = [cm_17_predict_noduplicates(:); cm_17_unpredict_noduplicates(:)];

num_subjects = 80; 
num_sessions1 = 7;
num_sessions2 = 7;
subj_num_mtx = repmat((1:num_subjects), num_sessions1*num_sessions2, 1); subj_num_mtx = subj_num_mtx(:); % subject matrix
session_num_mtx1 = repmat((1:num_sessions1)', num_sessions1, num_subjects); session_num_mtx1 = session_num_mtx1(:); % stage matrix
session_num_mtx2 = repmat(1:num_sessions2, num_sessions2, num_subjects); session_num_mtx2 = session_num_mtx2(:); % stage matrix
group_mtx = [zeros((num_subjects/2)*num_sessions1*num_sessions2, 1); ones((num_subjects/2)*num_sessions1*num_sessions2, 1)];

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
%var_name = 'ProbeBeh_adjacent_predict'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')

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
%var_name = 'ProbeBeh_adjacent_unpredict'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')

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
%var_name = 'ProbeBeh_adjacent_comp'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')

    %  stats: repeated measures anova for interaction
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
    
    
comp2seven_cell_unpredict









green_color = [21 154 72]./255; blue_color =  [43 73 159]./255;
light_green_color = [162 201 129]./255; light_blue_color = [102 164 213]./255;
figure_path = 'C:\Users\ampm1\Documents\manuscripts\gentrain\manuscript_figures\fig05\';
predict_subjs = [{'651049m1'} {'658648m2'} {'691359m1'} {'690330m1'} {'691359m2'} {'696944m3'} {'699437m4'} {'699438m3'}];
unpredict_subjs = [{'687034m4'} {'683472m2'} {'699437m2'} {'699438m1'} {'683472m3'}];

% for time warping
tses = [0.2 1.0 1.0 2.0 3.0 2.0];
event_frame = cumsum(tses(1:end-1)).*100; %from second np to reward delivery/ quit


%% example reliable cells

% all probe sessions
%{
[combine_trial_activity_mtx_predict, subject_ids_predict, session_nums_predict, cell_regist_mtx_cell_predict] = image_trial_activity_mtx_multisubj(predict_subjs, 1:3:19, 1:41, tses);
[combine_trial_activity_mtx_unpredict, subject_ids_unpredict, session_nums_unpredict, cell_regist_mtx_cell_unpredict] = image_trial_activity_mtx_multisubj(unpredict_subjs, 1:3:19, 1:41, tses);
save('combine_trial_activity_predict', 'combine_trial_activity_mtx_predict', 'subject_ids_predict', 'session_nums_predict', 'cell_regist_mtx_cell_predict')
save('combine_trial_activity_unpredict', 'combine_trial_activity_mtx_unpredict', 'subject_ids_unpredict', 'session_nums_unpredict', 'cell_regist_mtx_cell_unpredict')
%}

load('combine_trial_activity_predict'); load('combine_trial_activity_unpredict')

% sort index
[~,sort_idx_predict] = sort_rows_by_peak(norm_mtx(nanmean(combine_trial_activity_mtx_predict,3)')');
combine_trial_activity_mtx_predict_sorted = combine_trial_activity_mtx_predict(sort_idx_predict,:,:);

% iterate through cells of interest (see also 1240 2185)
for ic = [859 2126 2452]
    
    % plot
    figure; hold on
    imagesc(norm_mtx(permute(combine_trial_activity_mtx_predict_sorted(ic,:,:), [3 2 1])')'); 
    
    % aesthetics
    for ief = event_frame(1:end-1) 
       plot((ief-min(event_frame)+1).*[1 1], ylim, 'r-', 'linewidth', 1.0) 
    end
    xlim([0.5 size(combine_trial_activity_mtx_predict_sorted,2)+0.5])
    ylim([0.5 size(combine_trial_activity_mtx_predict_sorted,3)+0.5])
    set(gca, 'YDir','reverse'); xlabel('Event time'); ylabel('Trial')
    
    %save
    var_name = ['within_session_reliable_c' num2str(ic)]; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
end


%% distribution of firing fields over a trial
%{
[all_cell_corrs_predict, all_merge_mtx_predict, subj_corr_cell_predict, all_common_cell_matrices_predict, subj_corr_means_mtx_predict, subj_corr_cell_cellid_predict] = cell_turnover_timewarp_trials_multisubj('mevar', predict_subjs, 1:19, tses);
save('firing_field_data_predict', 'all_cell_corrs_predict', 'all_merge_mtx_predict', 'subj_corr_cell_predict', 'all_common_cell_matrices_predict', 'subj_corr_means_mtx_predict', 'subj_corr_cell_cellid_predict')
[all_cell_corrs_unpredict, all_merge_mtx_unpredict, subj_corr_cell_unpredict, all_common_cell_matrices_unpredict, subj_corr_means_mtx_unpredict, subj_corr_cell_cellid_unpredict] = cell_turnover_timewarp_trials_multisubj('hivar', unpredict_subjs, 1:19, tses);
save('firing_field_data_unpredict', 'all_cell_corrs_unpredict', 'all_merge_mtx_unpredict', 'subj_corr_cell_unpredict', 'all_common_cell_matrices_unpredict', 'subj_corr_means_mtx_unpredict', 'subj_corr_cell_cellid_unpredict')
%}
load('firing_field_data_predict')
load('firing_field_data_unpredict')

% single merged sorted matrix of all imaged cells (last appearance)
merge_mtx_predict = all_merge_mtx_predict(:,event_frame(1):event_frame(length(event_frame)-1));
    merge_mtx_predict = merge_mtx_predict(~isnan(sum(merge_mtx_predict,2)),:);
merge_mtx_unpredict = all_merge_mtx_unpredict(:,event_frame(1):event_frame(length(event_frame)-1));
    merge_mtx_unpredict = merge_mtx_unpredict(~isnan(sum(merge_mtx_unpredict,2)),:);
merge_mtx_all = [merge_mtx_predict; merge_mtx_unpredict];
merge_mtx_cell = [{merge_mtx_all}; {merge_mtx_predict}; {merge_mtx_unpredict}];
celltitles = {'All cells', 'Predictable training', 'Unpredictable training'};

% plot
figure;
for isp = 2:3
    
    % plot
    subplot(1,2,isp-1); hold on; 
    imagesc(sort_rows_by_peak(merge_mtx_cell{isp}))
    
    % aesthetics
    xticks_hold = [];
    for ief = event_frame 
       plot((ief-min(event_frame)+1).*[1 1], ylim, 'r-', 'linewidth', 1.0) 
       xticks_hold = [xticks_hold (ief-min(event_frame)+1)];
    end
    set(gca,'TickLength',[0, 0]); box off;
        title(celltitles{isp})
        ylabel('Neuron')
        xlabel('Event time')
        colorbar
    set(gca, 'YDir','reverse')
    ylim([0.5 size(merge_mtx_cell{isp},1)+0.5])
    xlim([0.5 size(merge_mtx_cell{isp},2)+0.5])
    xticks(xticks_hold)
    xticklabels_universe = {'NP', 'HE', 'TO', 'MW'};
    xticklabels(xticklabels_universe)
    yticks([1 size(merge_mtx_cell{isp},1)])
end

% save
set(gcf, 'Position', [6 399 1513 591])
var_name = 'cell_sorted_hms'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')



%% Held cells correlation matrices
load('firing_field_data_predict', 'subj_corr_means_mtx_predict')
load('firing_field_data_unpredict', 'subj_corr_means_mtx_unpredict')

figure;

subplot(1,2,1); 
imagesc(nanmean(subj_corr_means_mtx_predict,3)); 
xticks(1:3:19); xlabel('Imaging session'); yticks(1:3:19); ylabel('Imaging session');  
axis square; caxis([-.5 .75]); colorbar
set(gca,'TickLength',[0, 0]); title('Predictable')

subplot(1,2,2); 
imagesc(nanmean(subj_corr_means_mtx_unpredict,3)); 
xticks(1:3:19); xlabel('Imaging session'); yticks(1:3:19); ylabel('Imaging session');  
axis square; caxis([-.5 .75]); colorbar
set(gca,'TickLength',[0, 0]); title('Unpredictable')

var_name = 'held_cell_corr_hms'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')


save('Popcorr_raw', 'subj_corr_means_mtx_predict', 'subj_corr_means_mtx_unpredict')




%% pseudo session held cells
%{
[~, ~, ~, all_common_cell_pseudosesh_predict] = cell_turnover_timewarp_trials_multisubj('mevar', predict_subjs, [1 4 16 19], tses);
all_common_cell_pseudosesh_predict([2 4:end],:) = [];
[~, ~, ~, all_common_cell_pseudosesh_unpredict] = cell_turnover_timewarp_trials_multisubj('hivar', unpredict_subjs, [1 4 16 19], tses);
all_common_cell_pseudosesh_unpredict([2 4:end],:) = [];
save('psudosesh_figs', 'all_common_cell_pseudosesh_predict', 'all_common_cell_pseudosesh_unpredict'); 
%}
load('psudosesh_figs', 'all_common_cell_pseudosesh_predict', 'all_common_cell_pseudosesh_unpredict'); 

% plot
figure;

% predict probe 0 and 1
[~,sort_idx] = sort_rows_by_peak(all_common_cell_pseudosesh_predict{1,1}(:, event_frame(1):event_frame(length(event_frame)-1)));
subplot(2, 4, 1); imagesc(norm_mtx(all_common_cell_pseudosesh_predict{1,1}(sort_idx, event_frame(1):event_frame(length(event_frame)-1))')')
subplot(2, 4, 2); imagesc(norm_mtx(all_common_cell_pseudosesh_predict{1,2}(sort_idx, event_frame(1):event_frame(length(event_frame)-1))')')
% predict probe 5 and 6
[~,sort_idx] = sort_rows_by_peak(all_common_cell_pseudosesh_predict{2,1}(:, event_frame(1):event_frame(length(event_frame)-1)));
subplot(2, 4, 3); imagesc(norm_mtx(all_common_cell_pseudosesh_predict{2,1}(sort_idx, event_frame(1):event_frame(length(event_frame)-1))')')
subplot(2, 4, 4); imagesc(norm_mtx(all_common_cell_pseudosesh_predict{2,2}(sort_idx, event_frame(1):event_frame(length(event_frame)-1))')')
% unpredict probe 0 and 1
[~,sort_idx] = sort_rows_by_peak(all_common_cell_pseudosesh_unpredict{1,1}(:, event_frame(1):event_frame(length(event_frame)-1)));
subplot(2, 4, 5); imagesc(norm_mtx(all_common_cell_pseudosesh_unpredict{1,1}(sort_idx, event_frame(1):event_frame(length(event_frame)-1))')')
subplot(2, 4, 6); imagesc(norm_mtx(all_common_cell_pseudosesh_unpredict{1,2}(sort_idx, event_frame(1):event_frame(length(event_frame)-1))')')
% unpredict probe 5 and 6
[~,sort_idx] = sort_rows_by_peak(all_common_cell_pseudosesh_unpredict{2,1}(:, event_frame(1):event_frame(length(event_frame)-1)));
subplot(2, 4, 7); imagesc(norm_mtx(all_common_cell_pseudosesh_unpredict{2,1}(sort_idx, event_frame(1):event_frame(length(event_frame)-1))')')
subplot(2, 4, 8); imagesc(norm_mtx(all_common_cell_pseudosesh_unpredict{2,2}(sort_idx, event_frame(1):event_frame(length(event_frame)-1))')')

% aesthetics
for isp = 1:8
    subplot(2, 4, isp);hold on
    xticks_hold = [];
    for ief = event_frame(1:length(event_frame)-1)
       plot((ief-min(event_frame)+1).*[1 1], ylim, 'r-', 'linewidth', 1.0) 
       xticks_hold = [xticks_hold (ief-min(event_frame)+1)];
    end
    set(gca,'TickLength',[0, 0]); box off;
        ylabel('Neuron')
        xlabel('Event time')
    set(gca, 'YDir','reverse')
    ylim([0.5 inf])
    xlim([0.5 inf])
    xticks(xticks_hold)
    xticklabels_universe = {'NP', 'HE', 'TO', 'MW'};
    xticklabels(xticklabels_universe)
    yticks([1])
end

% save
set(gcf, 'Position', [47   294   990   684])
var_name = 'psudosesh_corr_eg_fig'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')



%% correct for time CORRELATION
load('days_between.mat', 'days_mtx_predict', 'days_mtx_unpredict')

min_cells = 1;

% average across cells in each session to compute subjects' mean
% correlation
subj_corr_means_mtx_predict = nan(19,19,length(predict_subjs)); 
subj_corr_ct_mtx_predict = nan(19,19,length(predict_subjs));
for isubj=1:length(predict_subjs)
    mtx_hold = subj_corr_means_mtx_predict(:,:,isubj);
    mtx_hold = cellfun(@mean, subj_corr_cell_predict{isubj});
    mtx_hold(logical(eye(length(mtx_hold)))) = nan;
    subj_corr_means_mtx_predict(:,:,isubj) = mtx_hold;
    subj_corr_ct_mtx_predict(:,:,isubj) = cellfun(@length, subj_corr_cell_predict{isubj});
end
subj_corr_means_mtx_unpredict = nan(19,19,length(unpredict_subjs)); 
subj_corr_ct_mtx_unpredict = nan(19,19,length(unpredict_subjs)); 
for isubj=1:length(unpredict_subjs)
    mtx_hold = subj_corr_means_mtx_unpredict(:,:,isubj);
    mtx_hold = cellfun(@mean, subj_corr_cell_unpredict{isubj});
    mtx_hold(logical(eye(length(mtx_hold)))) = nan;
    subj_corr_means_mtx_unpredict(:,:,isubj) = mtx_hold;
    subj_corr_ct_mtx_unpredict(:,:,isubj) = cellfun(@length, subj_corr_cell_unpredict{isubj});
end

% nan too few cells
subj_corr_means_mtx_predict(subj_corr_ct_mtx_predict < min_cells) = nan;
subj_corr_means_mtx_unpredict(subj_corr_ct_mtx_unpredict < min_cells) = nan;

% remove identity diaganol
%
subj_corr_means_mtx_predict(isinf(subj_corr_means_mtx_predict)) = nan; 
subj_corr_means_mtx_unpredict_smooth(isinf(subj_corr_means_mtx_unpredict)) = nan; 
subj_corr_means_mtx_predict(isinf(subj_corr_means_mtx_predict)) = nan; 
subj_corr_means_mtx_unpredict_smooth(isinf(subj_corr_means_mtx_unpredict)) = nan;
days_mtx_predict = days_mtx_predict(1:19,1:19,:);
    days_mtx_predict(days_mtx_predict==0) = nan;
days_mtx_unpredict = days_mtx_unpredict(1:19,1:19,:);
    days_mtx_unpredict(days_mtx_unpredict==0) = nan;
%}

% combine
%all_mean_activity_correlation = [subj_corr_means_mtx_predict(:); subj_corr_means_mtx_unpredict(:)]; 
    all_mean_activity_correlation = subj_corr_means_mtx_unpredict(:); 
%all_days = [days_mtx_predict(:);days_mtx_unpredict(:)];
    all_days = days_mtx_unpredict(:);
all_mean_activity_correlation = all_mean_activity_correlation(~isnan(all_days) & all_days<=90);
all_days = all_days(~isnan(all_days) & all_days<=90);
    
% model time
ibins = [0:0.25:5 5.5:.5:30 31:1:40 42:2:60 65:5:100];
ibin_cell = cell(1,length(ibins)-1);
for ibin = 1:length(ibin_cell)
    ibin_cell{ibin} = all_mean_activity_correlation(all_days>ibins(ibin) & all_days<=ibins(ibin+1));
end
max_x = 100;
[~, coefEsts, modelFun] = ampm_normal_fit_intzero(all_days(~isnan(all_days) & all_days<max_x), all_mean_activity_correlation(~isnan(all_days) & all_days<max_x), [1 0 150]);

%plot data points, boxcar mean, and fit curve
figure; hold on; 
plot(all_days, all_mean_activity_correlation, 'o', 'color', 0.8.*[1 1 1])
errorbar_plot_lineonly(ibin_cell, 0, ibins(1:end-1), 0.8.*[1 1 1])
plot(ibins, modelFun(coefEsts, ibins), 'r-')
%tic_vect = [-.999 -.99 -.9 -.8 -.6 -.3 0 .3 .6 .8 .9 .99 .999]; ylim(atanh([-.95 .95])); yticks(atanh(tic_vect)); yticklabels(tic_vect)
xlim([-1 90]); plot(xlim, [0 0], 'k--');
ylabel('Population activity correlation')
xlabel('Days between imaging sessions')
var_name = 'Time_curve_corr'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')

% correct population matrices
%
    % predict
    figure; 
    subplot(1,2,1); imagesc(nanmean(subj_corr_means_mtx_predict,3)); 
    axis square; caxis([-.5 .7]); colorbar; title('original'); set(gca,'TickLength',[0, 0]);
    xticks(1:2:19); xlabel('Imaging session'); yticks(1:2:19); ylabel('Imaging session')
    subj_corr_means_mtx_predict = subj_corr_means_mtx_predict - modelFun(coefEsts, days_mtx_predict);
    subplot(1,2,2); imagesc(nanmean(subj_corr_means_mtx_predict,3)); 
    axis square; caxis([-.4 .5]); colorbar; title('Time corrected'); set(gca,'TickLength',[0, 0]);
    xticks(1:2:19); xlabel('Imaging session'); yticks(1:2:19); ylabel('Imaging session')
    sgtitle Predict
    var_name = 'Corm_timecorrected_predict'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')

    % UNpredict
    figure; 
    subplot(1,2,1); imagesc(nanmean(subj_corr_means_mtx_unpredict,3)); 
    axis square; caxis([-.5 .7]); colorbar; title('original'); set(gca,'TickLength',[0, 0]);
    xticks(1:2:19); xlabel('Imaging session'); yticks(1:2:19); ylabel('Imaging session')
    subj_corr_means_mtx_unpredict = subj_corr_means_mtx_unpredict - modelFun(coefEsts, days_mtx_unpredict);
    subplot(1,2,2); imagesc(nanmean(subj_corr_means_mtx_unpredict,3)); 
    axis square; caxis([-.4 .5]); colorbar; title('Time corrected'); set(gca,'TickLength',[0, 0]);
    xticks(1:2:19); xlabel('Imaging session'); yticks(1:2:19); ylabel('Imaging session')
    sgtitle Unpredict
    var_name = 'Corm_timecorrected_unpredict'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')

    
% model comparing overlap matrices
subj_corr_means_mtx_predict_noduplicates = subj_corr_means_mtx_predict;
subj_corr_means_mtx_unpredict_noduplicates = subj_corr_means_mtx_unpredict;
for isesh1 = 1:size(subj_corr_means_mtx_predict,1)
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
lme_PreprobeFinalprobe = fitlme(tbl, model_str)

    model_str = 'MemorySimilarity~Session1*Session2+(1|Subject)';
    tbl = table(datamtx(group_mtx==0), session_num_mtx1(group_mtx==0), session_num_mtx2(group_mtx==0), subj_num_mtx(group_mtx==0), 'VariableNames', {'MemorySimilarity', 'Session1', 'Session2', 'Subject'});
    tbl.Subject = categorical(tbl.Subject);
    lme_PreprobeFinalprobe_predict = fitlme(tbl, model_str)
    
    model_str = 'MemorySimilarity~Session1*Session2+(1|Subject)';
    tbl = table(datamtx(group_mtx==1), session_num_mtx1(group_mtx==1), session_num_mtx2(group_mtx==1), subj_num_mtx(group_mtx==1), 'VariableNames', {'MemorySimilarity', 'Session1', 'Session2', 'Subject'});
    tbl.Subject = categorical(tbl.Subject);
    lme_PreprobeFinalprobe_unpredict = fitlme(tbl, model_str)


%% PopCorr comparisons (errorbars)
%{
min_cells = 2;
preprobe_session_numbers = 1:3:16;
firstday_session_numbers = 2:3:17;
lastday_session_numbers = 3:3:18;
postprobe_session_numbers = 4:3:19;

% preallcate subjects x probes
corr_predict_PreprobeProblem = nan(length(predict_subjs), length(preprobe_session_numbers));
corr_predict_ProblemFinalprobe = nan(length(predict_subjs), length(preprobe_session_numbers));
corr_unpredict_PreprobeProblem = nan(length(unpredict_subjs), length(preprobe_session_numbers));
corr_unpredict_ProblemFinalprobe = nan(length(unpredict_subjs), length(preprobe_session_numbers));
corr_predict_PreprobeFinalprobe = nan(length(predict_subjs), length(preprobe_session_numbers));
corr_unpredict_PreprobeFinalprobe = nan(length(unpredict_subjs), length(preprobe_session_numbers));
corr_predict_PreprobePostprobe = nan(length(predict_subjs), length(preprobe_session_numbers));
corr_unpredict_PreprobePostprobe = nan(length(unpredict_subjs), length(preprobe_session_numbers));

% iterate through sessions
for isesh = 1:length(preprobe_session_numbers)
    
    % preprobe and problem sessions
    %
        % iterate through predict subjects
        for isubj_predict = 1:length(predict_subjs)
            first_problem_rs = subj_corr_cell_predict{isubj_predict}{preprobe_session_numbers(isesh), firstday_session_numbers(isesh)};
                first_problem_rs = first_problem_rs - modelFun(coefEsts, days_mtx_predict(preprobe_session_numbers(isesh), firstday_session_numbers(isesh), isubj_predict));
            last_problem_rs = subj_corr_cell_predict{isubj_predict}{preprobe_session_numbers(isesh), lastday_session_numbers(isesh)};
                last_problem_rs = last_problem_rs - modelFun(coefEsts, days_mtx_predict(preprobe_session_numbers(isesh), lastday_session_numbers(isesh), isubj_predict));
            first_last_problem_rs = [first_problem_rs; last_problem_rs];
            if length(first_last_problem_rs) >= min_cells
                corr_predict_PreprobeProblem(isubj_predict, isesh) = mean(first_last_problem_rs);
            else
                corr_predict_PreprobeProblem(isubj_predict, isesh) = nan;
            end
        end
        % iterate through unpredict subjects
        for isubj_unpredict = 1:length(unpredict_subjs)
            first_problem_rs = subj_corr_cell_unpredict{isubj_unpredict}{preprobe_session_numbers(isesh), firstday_session_numbers(isesh)};
                first_problem_rs = first_problem_rs - modelFun(coefEsts, days_mtx_unpredict(preprobe_session_numbers(isesh), firstday_session_numbers(isesh), isubj_unpredict));
            last_problem_rs = subj_corr_cell_unpredict{isubj_unpredict}{preprobe_session_numbers(isesh), lastday_session_numbers(isesh)};
                last_problem_rs = last_problem_rs - modelFun(coefEsts, days_mtx_unpredict(preprobe_session_numbers(isesh), lastday_session_numbers(isesh), isubj_unpredict));
            first_last_problem_rs = [first_problem_rs; last_problem_rs];
            if length(first_last_problem_rs) >= min_cells
                corr_unpredict_PreprobeProblem(isubj_unpredict, isesh) = mean(first_last_problem_rs);
            else
                corr_unpredict_PreprobeProblem(isubj_unpredict, isesh) = nan;
            end
        end
        
    % compare all problems with last probe
    %
        % iterate through predict subjects
        for isubj_predict = 1:length(predict_subjs)
            first_problem_rs = subj_corr_cell_predict{isubj_predict}{postprobe_session_numbers(end), firstday_session_numbers(isesh)};
                first_problem_rs = first_problem_rs - modelFun(coefEsts, days_mtx_predict(postprobe_session_numbers(end), firstday_session_numbers(isesh), isubj_predict));
            last_problem_rs = subj_corr_cell_predict{isubj_predict}{postprobe_session_numbers(end), lastday_session_numbers(isesh)};
                last_problem_rs = last_problem_rs - modelFun(coefEsts, days_mtx_predict(postprobe_session_numbers(end), lastday_session_numbers(isesh), isubj_predict));
            first_last_problem_rs = [first_problem_rs; last_problem_rs]; 
            if length(first_last_problem_rs) >= min_cells
                corr_predict_ProblemFinalprobe(isubj_predict, isesh) = mean(first_last_problem_rs);
            else
                corr_predict_ProblemFinalprobe(isubj_predict, isesh) = nan;
            end
        end
        % iterate through unpredict subjects
        for isubj_unpredict = 1:length(unpredict_subjs)
            first_problem_rs = subj_corr_cell_unpredict{isubj_unpredict}{postprobe_session_numbers(end), firstday_session_numbers(isesh)};
                first_problem_rs = first_problem_rs - modelFun(coefEsts, days_mtx_unpredict(postprobe_session_numbers(end), firstday_session_numbers(isesh), isubj_unpredict));
            last_problem_rs = subj_corr_cell_unpredict{isubj_unpredict}{postprobe_session_numbers(end), lastday_session_numbers(isesh)};
                last_problem_rs = last_problem_rs - modelFun(coefEsts, days_mtx_unpredict(postprobe_session_numbers(end), lastday_session_numbers(isesh), isubj_unpredict));
            first_last_problem_rs = [first_problem_rs; last_problem_rs];
            if length(first_last_problem_rs) >= min_cells
                corr_unpredict_ProblemFinalprobe(isubj_unpredict, isesh) = mean(first_last_problem_rs);
            else
                corr_unpredict_ProblemFinalprobe(isubj_unpredict, isesh) = nan;
            end
        end
        
    % compare all probes with last probe
    %
        % iterate through predict subjects
        for isubj_predict = 1:length(predict_subjs)
            if subj_corr_ct_mtx_predict(preprobe_session_numbers(isesh), postprobe_session_numbers(isesh), isubj_predict) >= min_cells
                corr_predict_PreprobePostprobe(isubj_predict, isesh) = subj_corr_means_mtx_predict(preprobe_session_numbers(isesh), postprobe_session_numbers(isesh), isubj_predict);
            else
                corr_predict_PreprobePostprobe(isubj_predict, isesh) =  nan;
            end
        end
        % iterate through unpredict subjects
        for isubj_unpredict = 1:length(unpredict_subjs)
            if subj_corr_ct_mtx_unpredict(preprobe_session_numbers(isesh), postprobe_session_numbers(isesh), isubj_unpredict) >= min_cells
                corr_unpredict_PreprobePostprobe(isubj_unpredict, isesh) = subj_corr_means_mtx_unpredict(preprobe_session_numbers(isesh), postprobe_session_numbers(isesh), isubj_unpredict);
            else
                corr_unpredict_PreprobePostprobe(isubj_unpredict, isesh) = nan;
            end
        end
        
    % compare all preprobes with postprobes
    %
        % iterate through predict subjects
        for isubj_predict = 1:length(predict_subjs)
            if subj_corr_ct_mtx_predict(preprobe_session_numbers(isesh), postprobe_session_numbers(end), isubj_predict) >= min_cells
                corr_predict_PreprobeFinalprobe(isubj_predict, isesh) = subj_corr_means_mtx_predict(preprobe_session_numbers(isesh), postprobe_session_numbers(end), isubj_predict);
            else
                corr_predict_PreprobeFinalprobe(isubj_predict, isesh) =  nan;
            end
        end
        % iterate through unpredict subjects
        for isubj_unpredict = 1:length(unpredict_subjs)
            if subj_corr_ct_mtx_unpredict(preprobe_session_numbers(isesh), postprobe_session_numbers(end), isubj_unpredict) >= min_cells
                corr_unpredict_PreprobeFinalprobe(isubj_unpredict, isesh) = subj_corr_means_mtx_unpredict(preprobe_session_numbers(isesh), postprobe_session_numbers(end), isubj_unpredict);
            else
                corr_unpredict_PreprobeFinalprobe(isubj_unpredict, isesh) = nan;
            end
        end
end

% combine information for later mixed models
datamtx_PreprobeProblem = [corr_predict_PreprobeProblem; corr_unpredict_PreprobeProblem];
datamtx_PreprobeFinalprobe = [corr_predict_PreprobeFinalprobe; corr_unpredict_PreprobeFinalprobe];
datamtx_ProblemFinalprobe = [corr_predict_ProblemFinalprobe; corr_unpredict_ProblemFinalprobe];

% plot preprobe vs problem
%{
figure; hold on; 
errorbar_plot(num2cell(corr_predict_PreprobeProblem,1), 1, [], light_green_color, green_color)
errorbar_plot(num2cell(corr_unpredict_PreprobeProblem,1), 1, [], light_blue_color, blue_color)
xticks(1:size(corr_predict_PreprobeProblem,2)); xlabel('Probe number'); xticklabels(0:size(corr_predict_PreprobeFinalprobe,2)-1);
ylim([-.6 .7]); %ylim(atanh([-.7 .7])); yticks(atanh(-.9:.1:.9)); yticklabels(-.9:.1:.9)
ylabel('Mean activity correlation')
axis square; title('preprobe and problem')
hold on; plot(xlim, [0 0], 'k--');

    % mixed model
    num_subjects = length(predict_subjs) + length(unpredict_subjs); 
    num_sessions = length(preprobe_session_numbers);
    subj_num_mtx = repmat((1:num_subjects)', 1, num_sessions); % subject matrix
    session_num_mtx = repmat(1:num_sessions, num_subjects, 1); % stage matrix
    group_mtx = [zeros(length(predict_subjs), num_sessions); ones(length(unpredict_subjs), num_sessions)];
    model_str = 'PopCorr~Session*Group+(1|Subject)';
    tbl = table(datamtx_PreprobeProblem(:), session_num_mtx(:), group_mtx(:), subj_num_mtx(:), 'VariableNames', {'PopCorr', 'Session', 'Group', 'Subject'});
    tbl.Subject = categorical(tbl.Subject);
    lme_PreprobeProblem_errorbar = fitlme(tbl, model_str)

        % one way anova
        rmtbl_adj_PreprobeProblem_predict = simple_mixed_anova(corr_predict_PreprobeProblem);
        rmtbl_adj_PreprobeProblem_unpredict = simple_mixed_anova(corr_unpredict_PreprobeProblem);
        
    % ttest compare to zero
    pval_predict_PreprobeProblem_0 = nan(size(corr_predict_PreprobeProblem,2),1); 
        stats_predict_PreprobeProblem_0 = cell(size(corr_predict_PreprobeProblem,2),1);
    pval_unpredict_PreprobeProblem_0 = nan(size(corr_unpredict_PreprobeProblem,2),1); 
        stats_unpredict_PreprobeProblem_0 = cell(size(corr_unpredict_PreprobeProblem,2),1);
    for istage = 1:size(corr_predict_PreprobeProblem,2)
        [~,pval_predict_PreprobeProblem_0(istage),~,stats_predict_PreprobeProblem_0{istage}] = ttest(corr_predict_PreprobeProblem(:,istage));
        [~,pval_unpredict_PreprobeProblem_0(istage),~,stats_unpredict_PreprobeProblem_0{istage}] = ttest(corr_unpredict_PreprobeProblem(:,istage));
    end

    % ttest stats within groups
    pval_adj_predict_firstlast = nan(size(corr_predict_PreprobeProblem,2),1); stats_adj_predict_firstlast = cell(size(corr_predict_PreprobeProblem,2),1);
    pval_adj_unpredict_firstlast = nan(size(corr_unpredict_PreprobeProblem,2),1); stats_adj_unpredict_firstlast = cell(size(corr_unpredict_PreprobeProblem,2),1);
    for istage = 2:size(corr_predict_PreprobeProblem,2)
        [~,pval_adj_predict_firstlast(istage),~,stats_adj_predict_firstlast{istage}] = ttest(corr_predict_PreprobeProblem(:,1), corr_predict_PreprobeProblem(:,istage));
        [~,pval_adj_unpredict_firstlast(istage),~,stats_adj_unpredict_firstlast{istage}] = ttest(corr_unpredict_PreprobeProblem(:,1), corr_unpredict_PreprobeProblem(:,istage));
    end

    % ttest stats between groups
    pval_adj_comp_firstlast = nan(size(corr_unpredict_PreprobeProblem,2),1); stats_adj_comppval_adj_comp_firstlast = cell(size(corr_unpredict_PreprobeProblem,2),1);
    for istage = 1:size(corr_predict_PreprobeProblem,2)
        [~,pval_adj_comp_firstlast(istage),~,stats_adj_comppval_adj_comp_firstlast{istage}] = ttest2(corr_predict_PreprobeProblem(:,istage), corr_unpredict_PreprobeProblem(:,istage));
        pval_text(pval_adj_comp_firstlast(istage), istage-.4, 0.35)
    end
    
    var_name = 'PopCorr_PreprobeProblem'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')      
%}


% plot preprobe vs probe 6
%
figure; hold on; 
errorbar_plot(num2cell(corr_predict_PreprobeFinalprobe,1), 1, [], light_green_color, green_color)
errorbar_plot(num2cell(corr_unpredict_PreprobeFinalprobe,1), 1, [], light_blue_color, blue_color)
xticks(1:size(corr_predict_PreprobeFinalprobe,2)); xlabel('Probe number'); xticklabels(0:size(corr_predict_PreprobeFinalprobe,2)-1);
ylim([-.6 .7]); %ylim(atanh([-.7 .7])); yticks(atanh(-.9:.1:.9)); yticklabels(-.9:.1:.9)
ylabel('Mean activity correlation')
axis square; title('probes to last probe')
hold on; plot(xlim, [0 0], 'k--');


    % mixed model 
    num_subjects = length(predict_subjs) + length(unpredict_subjs); 
    num_sessions = length(preprobe_session_numbers);
    subj_num_mtx = repmat((1:num_subjects)', 1, num_sessions); % subject matrix
    session_num_mtx = repmat(1:num_sessions, num_subjects, 1); % stage matrix
    group_mtx = [zeros(length(predict_subjs), num_sessions); ones(length(unpredict_subjs), num_sessions)];
    model_str = 'PopCorr~Session*Group+(1|Subject)';
    tbl = table(datamtx_PreprobeFinalprobe(:), session_num_mtx(:), group_mtx(:), subj_num_mtx(:), 'VariableNames', {'PopCorr', 'Session', 'Group', 'Subject'});
    tbl.Subject = categorical(tbl.Subject);
    lme_PreprobeFinalprobe_errorbar = fitlme(tbl, model_str)
    
        % one way anova
        rmtbl_adj_PreprobeFinalprobe_predict = simple_mixed_anova(corr_predict_PreprobeFinalprobe);
        rmtbl_adj_PreprobeFinalprobe_unpredict = simple_mixed_anova(corr_unpredict_PreprobeFinalprobe);

    % ttest compare to zero
    pval_last_predict_0 = nan(size(corr_predict_PreprobeFinalprobe,2),1); stats_last_predict_0 = cell(size(corr_predict_PreprobeFinalprobe,2),1);
    pval_last_unpredict_0 = nan(size(corr_unpredict_PreprobeFinalprobe,2),1); stats_last_unpredict_0 = cell(size(corr_unpredict_PreprobeFinalprobe,2),1);
    for istage = 1:size(corr_predict_PreprobeFinalprobe,2)
        [~,pval_last_predict_0(istage),~,stats_last_predict_0{istage}] = ttest(corr_predict_PreprobeFinalprobe(:,istage));
        [~,pval_last_unpredict_0(istage),~,stats_last_unpredict_0{istage}] = ttest(corr_unpredict_PreprobeFinalprobe(:,istage));
    end

    % ttest stats within groups
    pval_last_predict = nan(size(corr_predict_PreprobeFinalprobe,2),1); stats_last_predict = cell(size(corr_predict_PreprobeFinalprobe,2),1);
    pval_last_unpredict = nan(size(corr_unpredict_PreprobeFinalprobe,2),1); stats_last_unpredict = cell(size(corr_unpredict_PreprobeFinalprobe,2),1);
    for istage = 2:size(corr_predict_PreprobeFinalprobe,2)
        [~,pval_last_predict(istage),~,stats_last_predict{istage}] = ttest(corr_predict_PreprobeFinalprobe(:,1), corr_predict_PreprobeFinalprobe(:,istage));
        [~,pval_last_unpredict(istage),~,stats_last_unpredict{istage}] = ttest(corr_unpredict_PreprobeFinalprobe(:,1), corr_unpredict_PreprobeFinalprobe(:,istage));
    end

    % ttest stats between group
    pval_last_comp = nan(size(corr_unpredict_PreprobeFinalprobe,2),1); stats_last_comp = cell(size(corr_unpredict_PreprobeFinalprobe,2),1);
    for istage = 1:size(corr_predict_PreprobeFinalprobe,2)
        [~,pval_last_comp(istage),~,stats_last_comp{istage}] = ttest2(corr_predict_PreprobeFinalprobe(:,istage), corr_unpredict_PreprobeFinalprobe(:,istage));
        pval_text(pval_last_comp(istage), istage-.4, 0.35)
    end
    
    %var_name = 'PopCorr_PreprobeFinalprobe'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
    
    
% plot problem vs probe 6
%
figure; hold on; 
errorbar_plot(num2cell(corr_predict_ProblemFinalprobe,1), 1, [], light_green_color, green_color)
errorbar_plot(num2cell(corr_unpredict_ProblemFinalprobe,1), 1, [], light_blue_color, blue_color)
xticks(1:size(corr_predict_ProblemFinalprobe,2)); xlabel('Problem number');
ylim([-.6 .7]); %ylim(atanh([-.7 .7])); yticks(atanh(-.9:.1:.9)); yticklabels(-.9:.1:.9)
ylabel('Mean activity correlation')
axis square; title('problem and last probe')
hold on; plot(xlim, [0 0], 'k--');


    % mixed model
    num_subjects = length(predict_subjs) + length(unpredict_subjs); 
    num_sessions = length(preprobe_session_numbers);
    subj_num_mtx = repmat((1:num_subjects)', 1, num_sessions); % subject matrix
    session_num_mtx = repmat(1:num_sessions, num_subjects, 1); % stage matrix
    group_mtx = [zeros(length(predict_subjs), num_sessions); ones(length(unpredict_subjs), num_sessions)];
    model_str = 'PopCorr~Session*Group+(1|Subject)';
    tbl = table(datamtx_ProblemFinalprobe(:), session_num_mtx(:), group_mtx(:), subj_num_mtx(:), 'VariableNames', {'PopCorr', 'Session', 'Group', 'Subject'});
    tbl.Subject = categorical(tbl.Subject);
    lme_ProblemFinalprobe_errorbar = fitlme(tbl, model_str)

        % one way anova
        rmtbl_adj_ProblemFinalprobe_predict = simple_mixed_anova(corr_predict_ProblemFinalprobe);
        rmtbl_adj_ProblemFinalprobe_unpredict = simple_mixed_anova(corr_unpredict_ProblemFinalprobe);

    % ttest compare to zero
    pval_last_predict_firstlast_0 = nan(size(corr_predict_ProblemFinalprobe,2),1); 
        stats_last_predict_firstlast_0 = cell(size(corr_predict_ProblemFinalprobe,2),1);
    pval_last_unpredict_firstlast_0 = nan(size(corr_unpredict_ProblemFinalprobe,2),1); 
        stats_unpredict_firstlast_0 = cell(size(corr_unpredict_ProblemFinalprobe,2),1);
    for istage = 1:size(corr_predict_ProblemFinalprobe,2)
        [~,pval_last_predict_firstlast_0(istage),~,stats_last_predict_firstlast_0{istage}] = ttest(corr_predict_ProblemFinalprobe(:,istage));
        [~,pval_last_unpredict_firstlast_0(istage),~,stats_unpredict_firstlast_0{istage}] = ttest(corr_unpredict_ProblemFinalprobe(:,istage));
    end

    % ttest stats within groups
    pval_last_predict_firstlast = nan(size(corr_predict_ProblemFinalprobe,2),1); 
        stats_last_predict_firstlast = cell(size(corr_predict_ProblemFinalprobe,2),1);
    pval_last_unpredict_firstlast = nan(size(corr_unpredict_ProblemFinalprobe,2),1); 
        stats_last_unpredict_firstlast = cell(size(corr_unpredict_ProblemFinalprobe,2),1);
    for istage = 2:size(corr_predict_ProblemFinalprobe,2)
        [~,pval_last_predict_firstlast(istage),~,stats_last_predict_firstlast{istage}] = ttest(corr_predict_ProblemFinalprobe(:,1), corr_predict_ProblemFinalprobe(:,istage));
        [~,pval_last_unpredict_firstlast(istage),~,stats_last_unpredict_firstlast{istage}] = ttest(corr_unpredict_ProblemFinalprobe(:,1), corr_unpredict_ProblemFinalprobe(:,istage));
    end

    % ttest stats between group
    pval_last_comp_firstlast = nan(size(corr_unpredict_ProblemFinalprobe,2),1); stats_last_comp_firstlast = cell(size(corr_unpredict_ProblemFinalprobe,2),1);
    for istage = 1:size(corr_predict_ProblemFinalprobe,2)
        [~,pval_last_comp_firstlast(istage),~,stats_last_comp_firstlast{istage}] = ttest2(corr_predict_ProblemFinalprobe(:,istage), corr_unpredict_ProblemFinalprobe(:,istage));
        pval_text(pval_last_comp_firstlast(istage), istage-.4, 0.35)
    end
    
    %var_name = 'PopCorr_ProblemFinalprobe'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')



%% Overlap comparisons (correlations)
%{
[~, ~, diff_waits_problem_preprobe_predict, diff_waits_problem_postprobe_predict, correlation_plot_cell_predict] = ALL_probe_vs_behavior_fit('train_mevar_imaging_hpc');
[~, ~, diff_waits_problem_preprobe_unpredict, diff_waits_problem_postprobe_unpredict, correlation_plot_cell_unpredict] = ALL_probe_vs_behavior_fit('train_hivar_imaging_hpc');
save('sfigure02_2', 'diff_waits_problem_preprobe_predict', 'diff_waits_problem_preprobe_unpredict', 'correlation_plot_cell_predict', 'correlation_plot_cell_unpredict')
%}
load('sfigure02_2', 'diff_waits_problem_preprobe_predict', 'diff_waits_problem_preprobe_unpredict', 'correlation_plot_cell_predict', 'correlation_plot_cell_unpredict')
%accuracy_of_prediction = [diff_waits_problem_preprobe_predict; diff_waits_problem_preprobe_unpredict];
cpcp_idx = 1;
first_day_discrimination = [correlation_plot_cell_predict{cpcp_idx}; correlation_plot_cell_unpredict{cpcp_idx}];
preprobe_problem_popcorr = [corr_predict_PreprobeProblem; corr_unpredict_PreprobeProblem];

%{
% First day discrimination vs Overlap (preprobe and problem (mean of first day and last day))  
model_str = 'PopCorr~Discrim+(1|Subject)';

    tbl = table(datamtx_PreprobeProblem(:), first_day_discrimination(:), subj_num_mtx(:), 'VariableNames', {'PopCorr', 'Discrim', 'Subject'});
    tbl.Subject = categorical(tbl.Subject);
    lme_PreprobeProblem_corr = fitlme(tbl, model_str);
   
    plot2var(correlation_plot_cell_predict{cpcp_idx}, corr_predict_PreprobeProblem, ...
                correlation_plot_cell_unpredict{cpcp_idx}, corr_unpredict_PreprobeProblem, 'First day discrimination and population overlap between preprobe and problem sessions')
            xlim([-5 5]); 
            ylim(atanh([-.7 .7])); yticks(atanh(-.9:.1:.9)); yticklabels(-.9:.1:.9)
            ylabel('Mean activity correlation')        
            xlabel('First day discrimination (d`)')
            hold on; plot(xlim, [0 0], 'k--'); plot([0 0], ylim, 'k--')
            var_name = 'FirstdayDiscrim_PopCorrPreprobeFirstLastproblem'; 
            print([figure_path var_name], '-dpdf', '-painters', '-bestfit')

            
% First day discrimination vs Overlap (preprobe and final probe)
model_str = 'PopCorr~Discrim+(1|Subject)';

    tbl = table(datamtx_ProblemFinalprobe(:), first_day_discrimination(:), subj_num_mtx(:), 'VariableNames', {'PopCorr', 'Discrim', 'Subject'});
    tbl.Subject = categorical(tbl.Subject);
    lme_PreprobeFinalprobe_corr = fitlme(tbl, model_str)

    plot2var(corr_predict_PreprobeProblem, corr_predict_PreprobeFinalprobe, ...
                corr_unpredict_PreprobeProblem, corr_unpredict_PreprobeFinalprobe, 'First day discrimination and population overlap between preprobe and last probe')
            xlim([-5 5]); 
            ylim(atanh([-.7 .7])); yticks(atanh(-.9:.1:.9)); yticklabels(-.9:.1:.9)
            ylabel('Mean activity correlation')
            xlabel('First day discrimination (d`)')
            hold on; plot(xlim, [0 0], 'k--'); plot([0 0], ylim, 'k--')
            var_name = 'FirstdayDiscrim_PopCorrPreprobeLastprobe'; 
            print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
            

% First day discrimination vs Overlap (problem and final probe)
model_str = 'PopCorr~Discrim+(1|Subject)';

    tbl = table(datamtx_PreprobeFinalprobe(:), first_day_discrimination(:), subj_num_mtx(:), 'VariableNames', {'PopCorr', 'Discrim', 'Subject'});
    tbl.Subject = categorical(tbl.Subject);
    lme_ProblemFinalprobe_corr = fitlme(tbl, model_str)

    plot2var(corr_predict_PreprobeProblem, corr_predict_ProblemFinalprobe, ...
                corr_unpredict_PreprobeProblem, corr_unpredict_ProblemFinalprobe, 'First day discrimination and population overlap between problem sessions and last probe')
            xlim([-5 5]); 
            ylim(atanh([-.7 .7])); yticks(atanh(-.9:.1:.9)); yticklabels(-.9:.1:.9)
            ylabel('Mean activity correlation')            
            ylabel('Population overlap')
            hold on; plot(xlim, [0 0], 'k--'); plot([0 0], ylim, 'k--')
            var_name = 'FirstdayDiscrim_PopCorrFirstLastproblemLastprobe'; 
            print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
%}

            
% Excess overlap with final probe vs excess overlap during training
%
% First day discrimination vs Overlap (preprobe and final probe)
model_str = 'PopCorr~Discrim+(1|Subject)';

    tbl = table(datamtx_PreprobeFinalprobe(:), preprobe_problem_popcorr(:), subj_num_mtx(:), 'VariableNames', {'PopCorr', 'Discrim', 'Subject'});
    tbl.Subject = categorical(tbl.Subject);
    lme_Reinstate_v_PreprobeFinalprobe_corr = fitlme(tbl, model_str)

    plot2var(corr_predict_PreprobeProblem, corr_predict_PreprobeFinalprobe, ...
                corr_unpredict_PreprobeProblem, corr_unpredict_PreprobeFinalprobe, 'Reinstatement and population corr between preprobe and last probe')
            axis([-1 1 -.6 .7]); xticks(-1:.5:1); yticks(-.6:.3:.6)
            ylabel('Mean activity correlation')
            xlabel('Reinstatement')
            hold on; plot(xlim, [0 0], 'k--'); plot([0 0], ylim, 'k--')
            var_name = 'Dotcorr_trainingCorr_PreprobeLastprobeCorr'; 
            %print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
            

% First day discrimination vs Overlap (problem and final probe)
model_str = 'PopCorr~Discrim+(1|Subject)';

    tbl = table(datamtx_ProblemFinalprobe(:), preprobe_problem_popcorr(:), subj_num_mtx(:), 'VariableNames', {'PopCorr', 'Discrim', 'Subject'});
    tbl.Subject = categorical(tbl.Subject);
    lme_Reinstate_v_ProblemFinalprobe_corr = fitlme(tbl, model_str)

    plot2var(corr_predict_PreprobeProblem, corr_predict_ProblemFinalprobe, ...
                corr_unpredict_PreprobeProblem, corr_unpredict_ProblemFinalprobe, 'Reinstatement and population corr between problem sessions and last probe')
            axis([-1 1 -.6 .7]); xticks(-1:.5:1); yticks(-.6:.3:.6)
            ylabel('Mean activity correlation')            
            ylabel('Reinstatement')
            hold on; plot(xlim, [0 0], 'k--'); plot([0 0], ylim, 'k--')
            var_name = 'Dotcorr_trainingCorr_ProblemLastprobeCorr'; 
            %print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
            
            

%% Internal functions
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
green_color = [21 154 72]./255; blue_color =  [43 73 159]./255;
light_green_color = [162 201 129]./255; light_blue_color = [102 164 213]./255;
figure_path = 'C:\Users\ampm1\Documents\manuscripts\gentrain\manuscript_figures\sfig04\';
predict_subjs = [{'651049m1'} {'658648m2'} {'691359m1'} {'690330m1'} {'691359m2'} {'696944m3'} {'699437m4'} {'699438m3'}];
unpredict_subjs = [{'687034m4'} {'683472m2'} {'699437m2'} {'699438m1'} {'683472m3'}];


%% cell origin color plot eg (footprints) 
%{
load('C:\Users\ampm1\Documents\MATLAB\generalization\cellreg_data\658648m2\cell_reg_658648m2.mat');
all_footprints_colored(cell_registered_struct.spatial_footprints_corrected, cell_registered_struct.cell_to_index_map)
var_name = 'Footprint_color_eg'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
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
%save('overlap_subj_matrices', 'predict_subjs', 'unpredict_subjs')

% plot individual subject probe overlap matrices
%
figure; 
for isubj = 1:length(predict_subjs)
subplot(1,length(predict_subjs),isubj);
imagesc(prop_predict(:,:,isubj)); caxis([0 .4]); axis square
end
sgtitle('predict')
var_name = 'Subjects_overlap_mtx_predict'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')

figure; 
for isubj = 1:length(unpredict_subjs)
subplot(1,length(unpredict_subjs),isubj);
imagesc(prop_unpredict(:,:,isubj)); caxis([0 .4]); axis square
end
sgtitle('unpredict')
var_name = 'Subjects_overlap_mtx_unpredict'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
%}
 

%% correct for time OVERLAP

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
load('days_between.mat', 'days_mtx_predict', 'days_mtx_unpredict')
load('overlap_subj_matrices', 'predict_subjs', 'unpredict_subjs')

% plot time between sessions
%
figure;
    % predict 
    subplot(1,2,1); imagesc(nanmean(days_mtx_predict,3))
    axis square; caxis([0 60]); colorbar; title('Predictable training')
    set(gca,'TickLength',[0, 0]);
    xticks(1:2:19); xlabel('Imaging session')
    yticks(1:2:19); ylabel('Imaging session')
    % unpredict
    subplot(1,2,2); imagesc(nanmean(days_mtx_unpredict,3))
    axis square; caxis([0 60]); colorbar; title('Unpredictable training')
    set(gca,'TickLength',[0, 0]);
    xticks(1:2:19); xlabel('Imaging session')
    yticks(1:2:19); ylabel('Imaging session')
    
    var_name = 'DaysBetween'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
    
% remove identity diaganol
prop_predict(prop_predict>0.99) = nan; 
prop_unpredict(prop_unpredict>0.99) = nan; 
days_mtx_predict = days_mtx_predict(1:19,1:19,:);
    days_mtx_predict(days_mtx_predict==0) = nan;
days_mtx_unpredict = days_mtx_unpredict(1:19,1:19,:);
    days_mtx_unpredict(days_mtx_unpredict==0) = nan;

% remove duplicates

    
% combine
max_x = 100;
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

%plot data points, boxcar mean, and fit curve
figure; hold on; 
plot(all_days, all_overlap, 'o', 'color', 0.8.*[1 1 1])
errorbar_plot_lineonly(ibin_cell, 0, ibins(1:end-1), 0.8.*[1 1 1])
plot(ibins, modelFun(coefEsts, ibins), 'r-')
ylim([-.01 1]); xlim([-1 90])
ylabel('Population overlap')
xlabel('Days between imaging sessions')
var_name = 'Time_curve_overlap'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')

% correct population overlap matrices
%
    % predict
    figure; 
    subplot(1,2,1); imagesc(nanmean(prop_predict,3)); 
    axis square; caxis([0 .4]); colorbar; title('original'); set(gca,'TickLength',[0, 0]);
    xticks(1:2:19); xlabel('Imaging session'); yticks(1:2:19); ylabel('Imaging session')
    prop_predict = prop_predict - modelFun(coefEsts, days_mtx_predict);
    subplot(1,2,2); imagesc(nanmean(prop_predict,3)); 
    axis square; caxis([-.15 .15]); colorbar; title('Time corrected'); set(gca,'TickLength',[0, 0]);
    xticks(1:2:19); xlabel('Imaging session'); yticks(1:2:19); ylabel('Imaging session')
    sgtitle('Predict overlap')
    var_name = 'Overlap_mtx_timecorrected_predict'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')

    % UNpredict
    figure; 
    subplot(1,2,1); imagesc(nanmean(prop_unpredict,3)); 
    axis square; caxis([0 .4]); colorbar; title('original'); set(gca,'TickLength',[0, 0]);
    xticks(1:2:19); xlabel('Imaging session'); yticks(1:2:19); ylabel('Imaging session')
    prop_unpredict = prop_unpredict - modelFun(coefEsts, days_mtx_unpredict);
    subplot(1,2,2); imagesc(nanmean(prop_unpredict,3)); 
    axis square; caxis([-.15 .15]); colorbar; title('Time corrected'); set(gca,'TickLength',[0, 0]);
    xticks(1:2:19); xlabel('Imaging session'); yticks(1:2:19); ylabel('Imaging session')
    sgtitle('Unpredict overlap')
    var_name = 'Overlap_mtx_timecorrected_unpredict'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
    
    
    
%% correct for time CORRELATION

load('firing_field_data_predict', 'subj_corr_means_mtx_predict')
load('firing_field_data_unpredict', 'subj_corr_means_mtx_unpredict')
load('days_between.mat', 'days_mtx_predict', 'days_mtx_unpredict')

% remove identity diaganol
subj_corr_means_mtx_predict(subj_corr_means_mtx_predict>0.99) = nan; 
subj_corr_means_mtx_unpredict(subj_corr_means_mtx_unpredict>0.99) = nan; 
days_mtx_predict = days_mtx_predict(1:19,1:19,:);
    days_mtx_predict(days_mtx_predict==0) = nan;
days_mtx_unpredict = days_mtx_unpredict(1:19,1:19,:);
    days_mtx_unpredict(days_mtx_unpredict==0) = nan;

% combine
%all_mean_activity_correlation = [subj_corr_means_mtx_predict(:); subj_corr_means_mtx_unpredict(:)]; 
    all_mean_activity_correlation = [subj_corr_means_mtx_unpredict(:)]; 
%all_days = [days_mtx_predict(:);days_mtx_unpredict(:)];
    all_days = [days_mtx_unpredict(:)];
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
ylim([-1 1]); xlim([-1 90]); plot(xlim, [0 0], 'k--');
ylabel('Population overlap')
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





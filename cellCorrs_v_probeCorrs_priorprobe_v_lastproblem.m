% cellCorrs_v_probeCorrs
% compares average correlation between shared cells with correlation
% between probe response profiles


% subject ids
subject_ids = [{'651049m1'} {'658648m2'} {'683472m2'} {'683472m3'}];
hivar_mevar = [{'mevar'} {'mevar'} {'hivar'} {'hivar'}];

%% compute correlation between all cells all sessions
%{
 [all_cell_corrs_mevar, all_merge_mtx_mevar, subj_corr_cell_mevar] = cell_turnover_timewarp_trials_multisubj('mevar', subject_ids([1 2]));
 [all_cell_corrs_hivar, all_merge_mtx_hivar, subj_corr_cell_hivar] = cell_turnover_timewarp_trials_multisubj('hivar', subject_ids([3 4]));
 save('cellCorrs_v_probeCorrs_prep.mat')
 %}
load('cellCorrs_v_probeCorrs_prep.mat')

% probe days only
prior_probe_days = [1:3:16];
last_problem_days = [3:3:18];
            
% minimum number of cells
min_cell = 5;

% average cell correlations
cell_corr_means_out_mevar = cell_corr_means(subj_corr_cell_mevar, prior_probe_days, last_problem_days, min_cell);
cell_corr_means_out_hivar = cell_corr_means(subj_corr_cell_hivar, prior_probe_days, last_problem_days, min_cell);

% matrices (session, subj)
cell_corr_mtx_mevar = nan(length(prior_probe_days), length(subj_corr_cell_mevar));
for isubj = 1:size(cell_corr_means_out_mevar,3)
    cell_corr_mtx_mevar(:,isubj) = diag(cell_corr_means_out_mevar(:,:,isubj));
end
cell_corr_mtx_hivar = nan(length(prior_probe_days), length(subj_corr_cell_hivar));
for isubj = 1:size(cell_corr_means_out_hivar,3)
    cell_corr_mtx_hivar(:,isubj) = diag(cell_corr_means_out_hivar(:,:,isubj));
end

cell_corr_mtx_mevar
cell_corr_mtx_hivar



%% compute population overlap all sessions
full_overlap_mtx_mevar = overlapping_population_average([{'651049m1'} {'658648m2'}]);
full_overlap_mtx_hivar = overlapping_population_average([{'683472m2'} {'683472m3'}]);

% probe days only
prior_probe_days = [1:3:16];
last_problem_days = [3:3:18];
            
% minimum number of cells
overlap_mtx_mevar = nan(length(prior_probe_days), size(full_overlap_mtx_mevar,3));
for isubj = 1:size(full_overlap_mtx_mevar,3)
    overlap_mtx_mevar(:,isubj) = diag(full_overlap_mtx_mevar(prior_probe_days, last_problem_days, isubj));
end

overlap_mtx_hivar = nan(length(prior_probe_days), size(full_overlap_mtx_hivar,3));
for isubj = 1:size(full_overlap_mtx_hivar,3)
    overlap_mtx_hivar(:,isubj) = diag(full_overlap_mtx_hivar(prior_probe_days, last_problem_days, isubj));
end

overlap_mtx_mevar
overlap_mtx_hivar



%% get preprobe and first problem wait times
[train_mean_waits_mevar, probe_model_waits_mevar] = probe_predict_lastproblem('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_imaging_hpc');
[train_mean_waits_hivar, probe_model_waits_hivar] = probe_predict_lastproblem('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_hivar_imaging_hpc');


% total error
%{
subj_mevar_error = nan(6,length(train_mean_waits_mevar));
for isubj_mevar = 1:length(train_mean_waits_mevar)
    subj_mevar_error(:,isubj_mevar) = sum(abs(train_mean_waits_mevar{isubj_mevar}-probe_model_waits_mevar{isubj_mevar}),2);
end
subj_hivar_error = nan(6,length(train_mean_waits_hivar));
for isubj_hivar = 1:length(train_mean_waits_hivar)
    subj_hivar_error(:,isubj_hivar) = sum(abs(train_mean_waits_hivar{isubj_hivar}-probe_model_waits_hivar{isubj_hivar}),2);
end

subj_mevar_error %(session, subj)
subj_hivar_error %(session, subj)
%}



% prediction strength (rich - poor)
%
% training r-p - preprobe r-p ; therefore, high values mean prediction
% improved, low values mean prediction got worse, extreme values mean
% prediction changed in general
%
subj_mevar_predstr = nan(6,length(train_mean_waits_mevar));
for isubj_mevar = 1:length(train_mean_waits_mevar)
    subj_mevar_predstr(:,isubj_mevar) = (train_mean_waits_mevar{isubj_mevar}(:,1)-train_mean_waits_mevar{isubj_mevar}(:,2)) - (probe_model_waits_mevar{isubj_mevar}(:,1)-probe_model_waits_mevar{isubj_mevar}(:,2));
end
subj_hivar_predstr = nan(6,length(train_mean_waits_hivar));
for isubj_hivar = 1:length(train_mean_waits_hivar)
    subj_hivar_predstr(:,isubj_hivar) = (train_mean_waits_hivar{isubj_hivar}(:,1)-train_mean_waits_hivar{isubj_hivar}(:,2)) - (probe_model_waits_hivar{isubj_hivar}(:,1)-probe_model_waits_hivar{isubj_hivar}(:,2));
end

% normalize
%subj_mevar_predstr = norm_mtx(subj_mevar_predstr);
%subj_hivar_predstr = norm_mtx(subj_hivar_predstr);

subj_mevar_predstr %(session, subj)
subj_hivar_predstr %(session, subj)



%% pre to post-probe coefficient changes
% get normal and log coefs for every subject and probe
[coefEsts_out_normal_mevar, coefEsts_out_log_mevar] = probe_coefs('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_imaging_hpc');
coefEsts_out_normal_mevar_delta = coefEsts_out_normal_mevar(2:end,:) - coefEsts_out_normal_mevar(1:end-1,:);
coefEsts_out_log_mevar_delta = coefEsts_out_log_mevar(2:end,:) - coefEsts_out_log_mevar(1:end-1,:);
[coefEsts_out_normal_hivar, coefEsts_out_log_hivar] = probe_coefs('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_hivar_imaging_hpc');
coefEsts_out_normal_hivar_delta = coefEsts_out_normal_hivar(2:end,:) - coefEsts_out_normal_hivar(1:end-1,:);
coefEsts_out_log_hivar_delta = coefEsts_out_log_hivar(2:end,:) - coefEsts_out_log_hivar(1:end-1,:);

% CORRECT LOG FOR DIRECTIONALITY
coefEsts_out_log_hivar_delta(1:2:end, :) = -coefEsts_out_log_hivar_delta(1:2:end, :);



%% plot correlation

% combine
combined_cellcorrs = [cell_corr_mtx_mevar(:);cell_corr_mtx_hivar(:)];
combined_prediction_strength_deltas = [subj_mevar_predstr(:);subj_hivar_predstr(:)];
combined_log_coefs = [coefEsts_out_log_mevar_delta(:);coefEsts_out_log_hivar_delta(:)];
combined_normal_coefs = [coefEsts_out_normal_mevar_delta(:);coefEsts_out_normal_hivar_delta(:)];

% bin cell corrs (atanh)
num_cell_corr_bins = 11; cell_corr_edges = linspace(-2,2,num_cell_corr_bins+1);
cell_corr_edges(1) = -inf; cell_corr_edges(end) = inf;
    [~,~,cell_corr_bins] = histcounts(combined_cellcorrs, cell_corr_edges);
    %[combined_cellcorrs cell_corr_bins]

% bin prediction deltas
num_pred_bins = 11; 
if min(combined_prediction_strength_deltas)==0 % if normalized
    pred_delta_edges = linspace(0,1,num_pred_bins+1);
    pred_delta_edges(1) = -inf; pred_delta_edges(end) = inf;
else
    pred_delta_edges = linspace(-15,15,num_pred_bins+1);
    pred_delta_edges(1) = -inf; pred_delta_edges(end) = inf;
end
    [~,~,pred_delta_bins] = histcounts(combined_prediction_strength_deltas, pred_delta_edges);
    %[combined_prediction_strength_deltas pred_delta_bins]


% normal matrix
plot_mtx_normal = nan(num_pred_bins, num_cell_corr_bins);

% iterate through prediction bins
for ipred_bin = 1:num_pred_bins

    % interate through cell corr bins
    for icell_corr_bin = 1:num_cell_corr_bins

        % load matrix
        plot_mtx_normal(ipred_bin, icell_corr_bin) = mean(combined_normal_coefs(pred_delta_bins==ipred_bin & cell_corr_bins==icell_corr_bin));
        
    end
end

figure; hold on; 
ampm_pcolor(plot_mtx_normal)
ylabel('Prediction strength change (s)')
yticks(1:1:num_pred_bins+2); yticklabels(round(pred_delta_edges*100)/100)
xlabel('Population state correlation (r)')
xticks(1:1:num_cell_corr_bins+2)
xticklabels((round(tanh(cell_corr_edges).*100)/100))
axis square
set(gca, 'YDir','normal')
axis([1 num_cell_corr_bins+1 1 num_pred_bins+1])
plot(((num_cell_corr_bins/2)+1).*[1 1], ylim, 'k--'); plot(xlim, ((num_cell_corr_bins/2)+1).*[1 1], 'k--')
colorbar; caxis([-10 10])
title('Normal delta')

% logistic matrix
plot_mtx_log = nan(num_pred_bins, num_cell_corr_bins);

% iterate through prediction bins
for ipred_bin = 1:num_pred_bins

    % interate through cell corr bins
    for icell_corr_bin = 1:num_cell_corr_bins

        % load matrix
        plot_mtx_log(ipred_bin, icell_corr_bin) = mean(combined_log_coefs(pred_delta_bins==ipred_bin & cell_corr_bins==icell_corr_bin));
        
    end
end

figure; hold on; 
ampm_pcolor(plot_mtx_log)
ylabel('Prediction strength change')
yticks(1:1:num_pred_bins+2); yticklabels(round(pred_delta_edges*100)/100)
xlabel('Population state correlation (r)')
xticks(1:1:num_cell_corr_bins+2)
xticklabels((round(tanh(cell_corr_edges).*100)/100))
axis square
set(gca, 'YDir','normal')
axis([1 num_cell_corr_bins+1 1 num_pred_bins+1])
plot(((num_cell_corr_bins/2)+1).*[1 1], ylim, 'k--'); plot(xlim, ((num_cell_corr_bins/2)+1).*[1 1], 'k--')
colorbar; caxis([-10 10])
title('Logistic delta')




%% plot overlap

% combine
combined_celloverlap = [overlap_mtx_mevar(:);overlap_mtx_hivar(:)];
combined_prediction_strength_deltas = [subj_mevar_predstr(:);subj_hivar_predstr(:)];
combined_log_coefs = [coefEsts_out_log_mevar_delta(:);coefEsts_out_log_hivar_delta(:)];
combined_normal_coefs = [coefEsts_out_normal_mevar_delta(:);coefEsts_out_normal_hivar_delta(:)];

% bin cell overlap
num_cell_overlap_bins = 11; cell_corr_edges = linspace(0, 0.5, num_cell_overlap_bins+1);
[~,~,cell_overlap_bins] = histcounts(combined_celloverlap, cell_corr_edges);
    %[combined_cellcorrs cell_corr_bins]

% bin prediction deltas
num_pred_bins = 11; 
if min(combined_prediction_strength_deltas)==0 % if normalized
    pred_delta_edges = linspace(0,1,num_pred_bins+1);
    pred_delta_edges(1) = -inf; pred_delta_edges(end) = inf;
else
    pred_delta_edges = linspace(-15,15,num_pred_bins+1);
    pred_delta_edges(1) = -inf; pred_delta_edges(end) = inf;
end
    [~,~,pred_delta_bins] = histcounts(combined_prediction_strength_deltas, pred_delta_edges);
    %[combined_prediction_strength_deltas pred_delta_bins]


% normal matrix
plot_mtx_normal = nan(num_pred_bins, num_cell_overlap_bins);

% iterate through prediction bins
for ipred_bin = 1:num_pred_bins

    % interate through cell corr bins
    for icell_overlap_bin = 1:num_cell_overlap_bins

        % load matrix
        plot_mtx_normal(ipred_bin, icell_overlap_bin) = mean(combined_normal_coefs(pred_delta_bins==ipred_bin & cell_overlap_bins==icell_overlap_bin));
        
    end
end

figure; hold on; 
ampm_pcolor(plot_mtx_normal)
ylabel('Prediction strength change (s)')
yticks(1:1:num_pred_bins+2); yticklabels(round(pred_delta_edges*100)/100)
xlabel('Population overlap (%)')
xticks(1:1:num_cell_overlap_bins+2)
xticklabels((round(tanh(cell_corr_edges).*100)/100))
axis square
set(gca, 'YDir','normal')
axis([1 num_cell_overlap_bins+1 1 num_pred_bins+1])
plot(((num_cell_overlap_bins/2)+1).*[1 1], ylim, 'k--'); plot(xlim, ((num_cell_overlap_bins/2)+1).*[1 1], 'k--')
colorbar; caxis([-10 10])
title('Normal delta')

% logistic matrix
plot_mtx_log = nan(num_pred_bins, num_cell_overlap_bins);

% iterate through prediction bins
for ipred_bin = 1:num_pred_bins

    % interate through cell corr bins
    for icell_overlap_bin = 1:num_cell_overlap_bins

        % load matrix
        plot_mtx_log(ipred_bin, icell_overlap_bin) = mean(combined_log_coefs(pred_delta_bins==ipred_bin & cell_overlap_bins==icell_overlap_bin));
        
    end
end

figure; hold on; 
ampm_pcolor(plot_mtx_log)
ylabel('Prediction strength change')
yticks(1:1:num_pred_bins+2); yticklabels(round(pred_delta_edges*100)/100)
xlabel('Population overlap (%)')
xticks(1:1:num_cell_overlap_bins+2)
xticklabels((round(tanh(cell_corr_edges).*100)/100))
axis square
set(gca, 'YDir','normal')
axis([1 num_cell_overlap_bins+1 1 num_pred_bins+1])
plot(((num_cell_overlap_bins/2)+1).*[1 1], ylim, 'k--'); plot(xlim, ((num_cell_overlap_bins/2)+1).*[1 1], 'k--')
colorbar; caxis([-10 10])
title('Logistic delta')


%% internal function

function cell_corr_means_out = cell_corr_means(subj_cell, sesh_nums_rows, sesh_nums_cols, min_cell)
% outputs a matrix subj_cell averages at the intersections defined by sesh_num inputs

    cell_corr_means_out = nan(length(sesh_nums_rows), length(sesh_nums_cols), length(subj_cell));
    for isubj = 1:length(subj_cell)
        for isesh1 = 1:length(sesh_nums_rows)
            for isesh2 = 1:length(sesh_nums_cols)
                if length(subj_cell{isubj}{isesh1, isesh2})>=min_cell
                    cell_corr_means_out(isesh1,isesh2,isubj) = nanmean(atanh(subj_cell{isubj}{sesh_nums_rows(isesh1), sesh_nums_cols(isesh1)}));
                else
                    cell_corr_means_out(isesh1,isesh2,isubj) = nan;
                end
            end
        end
    end

end

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
first_problem_days = [2:3:17];
sesh_idx = sort([prior_probe_days first_problem_days]);
%{
cell_corrs = [{subj_corr_cell_mevar{1}(sesh_idx(sesh_idx<=length(subj_corr_cell_mevar{1})), sesh_idx(sesh_idx<=length(subj_corr_cell_mevar{1})))}...
                {subj_corr_cell_mevar{2}(sesh_idx(sesh_idx<=length(subj_corr_cell_mevar{2})), sesh_idx(sesh_idx<=length(subj_corr_cell_mevar{2})))}...
                {subj_corr_cell_hivar{1}(sesh_idx(sesh_idx<=length(subj_corr_cell_hivar{1})), sesh_idx(sesh_idx<=length(subj_corr_cell_hivar{1})))}...
                {subj_corr_cell_hivar{2}(sesh_idx(sesh_idx<=length(subj_corr_cell_hivar{2})), sesh_idx(sesh_idx<=length(subj_corr_cell_hivar{2})))}];
%}
            
% minimum number of cells
min_cell = 1;

% average cell correlations
cell_corr_means_out_mevar = cell_corr_means(subj_corr_cell_mevar, prior_probe_days, first_problem_days, min_cell);
cell_corr_means_out_hivar = cell_corr_means(subj_corr_cell_hivar, prior_probe_days, first_problem_days, min_cell);

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

%% get preprobe and first problem wait times
[train_mean_waits_mevar, probe_model_waits_mevar] = probe_predict_firstproblem('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_imaging_hpc');
[train_mean_waits_hivar, probe_model_waits_hivar] = probe_predict_firstproblem('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_hivar_imaging_hpc');


% total error
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

subj_mevar_predstr %(session, subj)
subj_hivar_predstr %(session, subj)




%% plot

% figure
figure; hold on

% plot circles mevar
plot(cell_corr_mtx_mevar(:), subj_mevar_error(:), 'go')
plot(cell_corr_mtx_hivar(:), subj_hivar_error(:), 'bo')

% plot fit line
[r,p] = fit_line([cell_corr_mtx_mevar(:);cell_corr_mtx_hivar(:)], [subj_mevar_error(:);subj_hivar_error(:)], 0);
title(['r=' num2str(r) '; p=' num2str(round(p*10000)/10000)])
xlabel('Population firing similarity (r)')
ylabel('Memory inconsistency (s)') 

tic_vect = [-.9 -.8 -.6 -.3 0 .3 .6 .8 .9]; 
xlim(atanh([-.6 .92])); xticks(atanh(tic_vect)); xticklabels(tic_vect)
%ylim(atanh([-.6 .65])); yticks(atanh(tic_vect)); yticklabels(tic_vect)

hold on; plot([0 0], ylim, 'k--')
hold on; plot(xlim, [0 0], 'k--')

title('Total error')



% figure
figure; hold on

% plot circles mevar
plot(cell_corr_mtx_mevar(:), subj_mevar_predstr(:), 'go')
plot(cell_corr_mtx_hivar(:), subj_hivar_predstr(:), 'bo')

% plot fit line
[r,p] = fit_line([cell_corr_mtx_mevar(:);cell_corr_mtx_hivar(:)], [subj_mevar_predstr(:);subj_hivar_predstr(:)], 0);
title(['r=' num2str(r) '; p=' num2str(round(p*10000)/10000)])
xlabel('Population firing similarity (r)')
ylabel('Prediction error (s)') 

tic_vect = [-.9 -.8 -.6 -.3 0 .3 .6 .8 .9]; 
xlim(atanh([-.6 .92])); xticks(atanh(tic_vect)); xticklabels(tic_vect)
%ylim(atanh([-.6 .65])); yticks(atanh(tic_vect)); yticklabels(tic_vect)

hold on; plot([0 0], ylim, 'k--')
hold on; plot(xlim, [0 0], 'k--')

title('Prediction strength error')


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

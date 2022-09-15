
% PROBE 8

green_color = [21 154 72]./255; blue_color =  [43 73 159]./255;
light_green_color = [162 201 129]./255; light_blue_color = [102 164 213]./255;
figure_path = 'C:\Users\ampm1\Documents\manuscripts\gentrain\manuscript_figures\fig05\';
predict_subjs = [{'651049m1'} {'658648m2'} {'691359m1'} {'690330m1'} {'691359m2'} {'696944m3'} {'699437m4'} {'699438m3'}];
unpredict_subjs = [{'687034m4'} {'683472m2'} {'699437m2'} {'699438m1'} {'683472m3'}];

% for time warping
tses = [0.2 1.0 1.0 2.0 3.0 2.0];
event_frame = cumsum(tses(1:end-1)).*100; %from second np to reward delivery/ quit

%% unoperated

%% Ca probe 6 and 7 similarity
figure;

% predict
probe_fps_predict = [get_file_paths_targeted('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_imaging_hpc', 'probe')];            
probe_paths_predict = cell(1, 2);
for iprobe = 1:2
    
    if iprobe==1
        probe_paths_predict{iprobe} = probe_fps_predict(contains(probe_fps_predict, 'postprobe_01'));
    elseif iprobe==2
        probe_paths_predict{iprobe} = probe_fps_predict(contains(probe_fps_predict, 'postprobe_02'));
    end
    
    subplot(2,2,iprobe)
    plot_allprobes(probe_paths_predict{iprobe}); 
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
    title(['Probe ' num2str(iprobe+6)])
end

% unpredict
probe_fps_unpredict = [get_file_paths_targeted('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_hivar_imaging_hpc', 'probe')];            
probe_paths_unpredict = cell(1, 2);
for iprobe = 1:2
    
    if iprobe==1
        probe_paths_unpredict{iprobe} = probe_fps_unpredict(contains(probe_fps_unpredict, 'postprobe_01'));
    elseif iprobe==2
        probe_paths_unpredict{iprobe} = probe_fps_unpredict(contains(probe_fps_unpredict, 'postprobe_02'));
    end
    
    subplot(2,2,iprobe+2)
    plot_allprobes(probe_paths_unpredict{iprobe}); 
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
end



%% probe correlations

interp_num = 4;

% predict
probe_corrs_predict = nan(length(predict_subjs), 1);
for isubj = 1:size(probe_paths_predict{1}, 1)
    
    % probe 6
    load(probe_paths_predict{1}{isubj}, 'trl_mtx')
    [wait_times_p6] = wait_times_prep(trl_mtx, 2, interp_num);
    
    % probe 7
    load(probe_paths_predict{2}{isubj}, 'trl_mtx')
    [wait_times_p7] = wait_times_prep(trl_mtx, 2, interp_num);
    
    if size(wait_times_p6,2)>size(wait_times_p6,1)
        wait_times_p6 = wait_times_p6'; wait_times_p7 = wait_times_p7';
    end
    nnan_idx = ~isnan(wait_times_p6) & ~isnan(wait_times_p7);
    probe_corrs_predict(isubj) = corr(wait_times_p6(nnan_idx), wait_times_p7(nnan_idx));
    
    figure; hold on; 
    plot(wait_times_p6(nnan_idx), 'ro')
    plot(wait_times_p7(nnan_idx), 'bo')
    title([find_subj_id(probe_paths_predict{1}{isubj}), '; r = ' num2str(probe_corrs_predict(isubj))])

end

% unpredict
probe_corrs_unpredict = nan(length(unpredict_subjs), 1);
for isubj = 1:size(probe_paths_unpredict{1}, 1)
    
    % probe 6
    load(probe_paths_unpredict{1}{isubj}, 'trl_mtx')
    [wait_times_p6] = wait_times_prep(trl_mtx, 2, interp_num);
    
    % probe 7
    load(probe_paths_unpredict{2}{isubj}, 'trl_mtx')
    [wait_times_p7] = wait_times_prep(trl_mtx, 2, interp_num);
    
    if size(wait_times_p6,2)>size(wait_times_p6,1)
        wait_times_p6 = wait_times_p6'; wait_times_p7 = wait_times_p7';
    end
    nnan_idx = ~isnan(wait_times_p6) & ~isnan(wait_times_p7);
    probe_corrs_unpredict(isubj) = corr(wait_times_p6(nnan_idx), wait_times_p7(nnan_idx));
    
    figure; hold on; 
    plot(wait_times_p6(nnan_idx), 'ro')
    plot(wait_times_p7(nnan_idx), 'bo')
    title([find_subj_id(probe_paths_unpredict{1}{isubj}), '; r = ' num2str(probe_corrs_unpredict(isubj))])
    
end

figure
errorbar_barplot([{probe_corrs_predict} {probe_corrs_unpredict}],[], [], [green_color; blue_color])
[~, p, ~, stats] = ttest2(probe_corrs_predict, probe_corrs_unpredict);
ylabel('Wait times similarity (r)')
ylim([-1 1])
title(['t(' num2str(stats.df) ') = ' num2str(round(stats.tstat*100)/100) '; p = ' num2str(round(p*100)/100)])
xticklabels({'Predict', 'Unpredict'})


%% population overlap

% compute cell registration matrices
%{
session_chron_reorder_idx = [1 9 15 2 10 16 3 11 17 4 12 18 5 13 19 6 14 20 7 8];
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
save('crm_1920', 'subject_cell_crm_predict', 'subject_cell_crm_unpredict')
%}
load('crm_1920', 'subject_cell_crm_predict', 'subject_cell_crm_unpredict')

session_numbers = 19:20;

% preallcate sessions x sessions x subjects
popoverlap_predict = nan(length(predict_subjs), 1);
popoverlap_unpredict = nan(length(unpredict_subjs), 1);

% iterate through predict subjects
for isubj_predict = 1:length(predict_subjs)
    % load crm
    crm = subject_cell_crm_predict{isubj_predict};
    % number unique cells
    unq_ct = sum(~isnan(crm(:,session_numbers(1))) | ~isnan(crm(:,session_numbers(2))));
    % active in both
    common_ct = sum(~isnan(crm(:,session_numbers(1))) & ~isnan(crm(:,session_numbers(2))));
    % load 
    popoverlap_predict(isubj_predict) = common_ct/unq_ct;
end

% iterate through unpredict subjects
for isubj_unpredict = 1:length(unpredict_subjs)
    % load crm
    crm = subject_cell_crm_unpredict{isubj_unpredict};
    % number unique cells
    unq_ct = sum(~isnan(crm(:,session_numbers(1))) | ~isnan(crm(:,session_numbers(2))));
    % active in both
    common_ct = sum(~isnan(crm(:,session_numbers(1))) & ~isnan(crm(:,session_numbers(2))));
    % load 
    popoverlap_unpredict(isubj_unpredict) = common_ct/unq_ct;
end

%figure
%errorbar_barplot([{popoverlap_predict} {popoverlap_unpredict}])
%[~, p, ~, stats] = ttest2(popoverlap_predict, popoverlap_unpredict);


%figure; hold on;
%plot(popoverlap_predict, probe_corrs_predict, 'go')
%plot(popoverlap_unpredict, probe_corrs_unpredict, 'bo')
%[r, p] = fit_line([popoverlap_predict; popoverlap_unpredict], [probe_corrs_predict; probe_corrs_unpredict], 0)



%% activity correlation
%{
[all_cell_corrs_predict, all_merge_mtx_predict, subj_corr_cell_predict, all_common_cell_matrices_predict, subj_corr_means_mtx_predict, subj_corr_cell_cellid_predict] = cell_turnover_timewarp_trials_multisubj('mevar', predict_subjs, 19:20, tses);
save('firing_field_data_predict_1920', 'all_cell_corrs_predict', 'all_merge_mtx_predict', 'subj_corr_cell_predict', 'all_common_cell_matrices_predict', 'subj_corr_means_mtx_predict', 'subj_corr_cell_cellid_predict')
[all_cell_corrs_unpredict, all_merge_mtx_unpredict, subj_corr_cell_unpredict, all_common_cell_matrices_unpredict, subj_corr_means_mtx_unpredict, subj_corr_cell_cellid_unpredict] = cell_turnover_timewarp_trials_multisubj('hivar', unpredict_subjs, 19:20, tses);
save('firing_field_data_unpredict_1920', 'all_cell_corrs_unpredict', 'all_merge_mtx_unpredict', 'subj_corr_cell_unpredict', 'all_common_cell_matrices_unpredict', 'subj_corr_means_mtx_unpredict', 'subj_corr_cell_cellid_unpredict')
%}
load('firing_field_data_predict_1920', 'subj_corr_cell_predict')
load('firing_field_data_unpredict_1920', 'subj_corr_cell_unpredict')

popcorr_predict = nan(length(predict_subjs), 1);
for isubj_predict = 1:length(subj_corr_cell_predict)
   popcorr_predict(isubj_predict) = nanmean(subj_corr_cell_predict{isubj_predict}{1,2}) ;
end
popcorr_unpredict = nan(length(unpredict_subjs), 1);
for isubj_unpredict = 1:length(subj_corr_cell_unpredict)
   popcorr_unpredict(isubj_unpredict) = nanmean(subj_corr_cell_unpredict{isubj_unpredict}{1,2}) ;
end

%figure
%errorbar_barplot([{popcorr_predict} {popcorr_unpredict}])
%[~, p, ~, stats] = ttest2(popcorr_predict, popcorr_unpredict);


%figure; hold on;
%plot(popcorr_predict, probe_corrs_predict, 'go')
%plot(popcorr_unpredict, probe_corrs_unpredict, 'bo')
%[r, p] = fit_line([popcorr_predict; popcorr_unpredict], [probe_corrs_predict; probe_corrs_unpredict], 0)


%% reinstatement score
r_score_predict = popoverlap_predict + popoverlap_predict .* popcorr_predict;
r_score_unpredict = popoverlap_unpredict + popoverlap_unpredict .* popcorr_unpredict;

figure
errorbar_barplot([{r_score_predict} {r_score_unpredict}], [], [], [green_color;blue_color])
title('Reactivation')
[~, p, ~, stats] = ttest2(r_score_predict, r_score_unpredict);
ylabel('Reactivation score')
ylim([0 1])
title(['t(' num2str(stats.df) ') = ' num2str(round(stats.tstat*100)/100) '; p = ' num2str(round(p*100)/100)])
xticklabels({'Predict', 'Unpredict'})

figure; hold on;
title('Reactivation')
plot(r_score_predict, probe_corrs_predict, 'o', 'color', green_color)
plot(r_score_unpredict, probe_corrs_unpredict, 'o', 'color', blue_color)
[r, p] = fit_line([r_score_predict; r_score_unpredict], [probe_corrs_predict; probe_corrs_unpredict], 0);
xlabel('Reactivation score')
ylabel('Wait times similarity (r)')
title(['r = ' num2str(round(r*100)/100) '; p = ' num2str(round(p*100)/100)])



%% internal functions

% sort combined crms into one psuedo subject matrix
function crm = treat_crm(crm, session_chron_reorder_idx)

    crm = crm(:, session_chron_reorder_idx);
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
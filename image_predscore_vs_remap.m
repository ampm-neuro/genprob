function [accuracy_of_prediction, pop_overlap, remap_score] = image_predscore_vs_remap(training_group)
% Compute prediction score of prior probe (to subsequent problem)
% and then quantify hpc remap from that probe to the next


%% unique subject folders in training folder

% could be input
folderpath = ['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\' training_group];

% file list
file_list_subjects = dir(folderpath);
file_list_subjects(1:2) = [];

% unique subject strings
unq_subjs = cell(1, length(file_list_subjects));
for isubj = 1:length(file_list_subjects)
    unq_subjs{isubj} = file_list_subjects(isubj).name;
end



%% Problem tones for this training group

% all tones
load('unqfrq41', 'unqfrq41')

% tones in each problem (rich, poor)
if contains(training_group, 'mevar')
    all_prob_tones = rich_bounds_prob('mevar', 0);
elseif contains(training_group, 'hivar')
    all_prob_tones = rich_bounds_prob('hivar', 0);
else
    error('input error')
end

% number of problems
num_problems = size(all_prob_tones,1);

% tones numbers (1:41) in each problem (rich, poor)
all_prob_tone_nums = nan(size(all_prob_tones));
for iprob = 1:num_problems
    for itone = 1:size(all_prob_tones,2)
        all_prob_tone_nums(iprob, itone) = find(ismember(unqfrq41, all_prob_tones(iprob, itone)));
    end
end



%% Discrimination problem behavior

% all wait times over learning
[firstLast_dprime_mtx, firstLast_waitTimes_cell, days_to_crit, all_TwoDay_learn_mtx] = ALL_stage_learn(['two_tone\' training_group], 1:6, 2);
%close; close; close; close; close; close; close
% firstLast_waitTimes_cell:
%   6 cells (1 per problem)
%       2 cells per subj (rows: subj, columns: firstDay, lastDay)
%           2 cells (rich, poor)
%               vector of wait times

% first and last day dprimes (subj, problem)
first_day_d = squeeze(firstLast_dprime_mtx(:,1,:));
last_day_d = squeeze(firstLast_dprime_mtx(:,2,:));

% total learning (last day performance minus first day performance)
lastMinusFirst_day_d = last_day_d-first_day_d;

% learning over first two days
TwoDay_learn_mtx = nan(size(all_TwoDay_learn_mtx,1), num_problems);
for iprob = 1:num_problems
    TwoDay_learn_mtx(:,iprob) = all_TwoDay_learn_mtx(:, 2, iprob)-all_TwoDay_learn_mtx(:, 1, iprob);
end

%% all probe wait times

% smooth wait times from every probe for every subject (probe,tone,subj)
%[probe_wait_times, unq_subjects_probe_all, subj_path_cell] = ALL_probe_wait_times(training_group, 0);
[probe_wait_times, unq_subjects_probe_all, subj_path_cell] = ALL_probe_wait_times(training_group, 1);
%[probe_wait_times, unq_subjects_probe_all, subj_path_cell] = ALL_probe_wait_times(training_group, 2);

% number of probes
num_probes = size(probe_wait_times,1);
num_tones = size(probe_wait_times,2);
num_subjects = size(probe_wait_times,3);

% zscored probe wait times
probe_wait_times_z = nan(num_probes, num_tones, num_subjects);
for iprobe = 1:num_probes
    for isubj = 1:num_subjects
        probe_wait_times_z(iprobe, :, isubj) = zscore_mtx(probe_wait_times(iprobe, :, isubj)')';
    end
end



%% mixed model prep
subj_num_mtx = repmat((1:num_subjects)', 1, num_problems); % subject matrix
stage_num_mtx = repmat(1:6, num_subjects, 1); % stage matrix
model_str = 'var1~var2+(1|subject)+(1|problem)';
model_str_against = 'var1~var2+problem+(1|subject)';
model_str_against_double = 'var1~var2+var3+problem+(1|subject)';



%% Fit model to each probe session (wait times)
%
% preallcoate 
num_coefs = 4;
coefficients = nan(num_probes, num_coefs, num_subjects);
subj_curves = nan(num_probes, length(1:0.001:num_tones), num_subjects);

% iterate through all probes
for iprobe = 1:num_probes
    
    % open curve of all subject curves for this probe
    figure; hold on

    % iterate through subjects
    for isubj = 1:num_subjects
        
        % isolate probe waits
        wait_times = probe_wait_times(iprobe, :, isubj);

        
        isubj
        [1:num_tones;
            wait_times]
        [nanmean(probe_wait_times(iprobe,:,isubj)) 20 0 1]
        
        % compute coefficients
        try
            mean_wait_times = nanmean(probe_wait_times(iprobe,:,isubj));
            %[~, coefEsts, modelFun] = ampm_normal_logistic_fit(1:num_tones, wait_times, [mean_wait_times 21 0 1]);
            [~, coefEsts, modelFun] = ampm_normal_logistic_fit_algo(1:num_tones, wait_times, [mean_wait_times 20 0 1]);
            
        % if modelling fails, set to mean
        catch 
            
            coefEsts = [nanmean(probe_wait_times(iprobe,:,isubj)) zeros(1,num_coefs-1)];
        end
        
        % load coefficients
        coefficients(iprobe,:,isubj) = coefEsts;
        
        % santity check coefs
        %acceptable_coef_bounds = [-5 45; -5 41; -20 20; -20 20];
        acceptable_coef_bounds = [-30 60; 1 41; -30 30; -30 30];
        if any(...
                coefficients(iprobe,1,isubj) < acceptable_coef_bounds(1,1) | coefficients(iprobe,1,isubj) > acceptable_coef_bounds(1,2) | isnan(coefficients(iprobe,1,isubj)) | ...
                coefficients(iprobe,2,isubj) < acceptable_coef_bounds(2,1) | coefficients(iprobe,2,isubj) > acceptable_coef_bounds(2,2) | isnan(coefficients(iprobe,2,isubj)) |...
                coefficients(iprobe,3,isubj) < acceptable_coef_bounds(3,1) | coefficients(iprobe,3,isubj) > acceptable_coef_bounds(3,2) | isnan(coefficients(iprobe,3,isubj)) |...
                coefficients(iprobe,4,isubj) < acceptable_coef_bounds(4,1) | coefficients(iprobe,4,isubj) > acceptable_coef_bounds(4,2) | isnan(coefficients(iprobe,4,isubj)) ...
                )
            coefficients(iprobe,:,isubj) = [nanmean(probe_wait_times(iprobe,:,isubj)) zeros(1,num_coefs-1)];
        end
        
        % deal with missing data
        if isnan(coefficients(iprobe,1,isubj))
            coefficients(iprobe,:,isubj) = nan;
        end
        
            % plot subject curve
            %
            iprobe
            %if iprobe == 7
            figure; hold on
            %wait_times
            plot(1:num_tones, wait_times,'o')
            subj_curves(iprobe,:,isubj) = modelFun(coefEsts, (1:0.001:num_tones)');
            plot(1:0.001:num_tones, subj_curves(iprobe,:,isubj), '-', 'color', 0.7.*[1 1 1], 'linewidth', 1)
            drawnow
            %coefficients(iprobe,:,isubj)
            %error
            title(['subj ' num2str(isubj) '; probe ' num2str(iprobe)])
            %end
            %}
    end
    
    % plot mean and se of all subject curves
    figure
    mean_fit_curve = nanmean(subj_curves(iprobe,:,:),3);
    se_fit_curves = nanstd(subj_curves(iprobe,:,:),[],3)./sqrt(sum(~isnan(subj_curves(iprobe,1,:)),3));
    plot(1:0.001:num_tones, mean_fit_curve, 'k-', 'linewidth', 4)
    plot(1:0.001:num_tones, mean_fit_curve-se_fit_curves, 'k-', 'linewidth', 2)
    plot(1:0.001:num_tones, mean_fit_curve+se_fit_curves, 'k-', 'linewidth', 2)
    set(gca,'TickLength',[0, 0]); box off;
    xlim([0.5 41.5])
    ylim([0 40])
    ylabel('Wait times (s)')
    xlabel('Tone number')
    title(['Probe ' num2str(iprobe) ': Norm = ' num2str(mean(coefficients(iprobe,4,:))) ', Logistic = ' num2str(mean(coefficients(iprobe,3,:)))])
    %}
    
end



%% Plot model coefficients

% Logistic curve mulitplier
logist_ceof = squeeze(coefficients(:,3,:))';
logist_ceof_corrected = squeeze(coefficients(:,3,:))';
logist_coef_cell = cell(1,num_probes);
logist_coef_cell_corrected = cell(1,num_probes);
for iprobe = 1:num_probes
    if ismember(iprobe, [2 4 6])
        logist_coef_cell_corrected{iprobe} = -logist_ceof(:,iprobe);
        logist_ceof_corrected(:,iprobe) = -logist_ceof(:,iprobe);
    else
        logist_coef_cell_corrected{iprobe} = logist_ceof(:,iprobe);
    end
    logist_coef_cell{iprobe} = logist_ceof(:,iprobe);
end
%{
figure; hold on;
errorbar_plot(logist_coef_cell_corrected(1), 0, 1);
errorbar_plot(logist_coef_cell_corrected(2:7), 1, 2:7);
errorbar_plot(logist_coef_cell_corrected(8), 0, 8);
errorbar_plot(logist_coef_cell_corrected(9), 0, 9);
%errorbar_plot(logist_coef_cell(1), 0, 1);
%errorbar_plot(logist_coef_cell(2:7), 1, 2:7);
%errorbar_plot(logist_coef_cell(8), 0, 8);
%errorbar_plot(logist_coef_cell(9), 0, 9);

%tests
for i = 1:9
    disp(['probe ' num2str(i)])
   [~,coef_log_p] =  ttest(logist_coef_cell{i});
end

xlim([0.5 num_probes+0.5])
plot(xlim, [1 1].*0, 'k--')
set(gca,'TickLength',[0, 0]); box off;
xlabel('Probe number')
ylabel('Wait times (s)')
title('Logistic curve height')
%}

% Normal curve mulitplier
norm_coef = squeeze(coefficients(:,4,:))';
norm_coef_cell = cell(1,num_probes);
for iprobe = 1:num_probes
    norm_coef_cell{iprobe} = norm_coef(:,iprobe);
end
%{
figure; hold on;
errorbar_plot(norm_coef_cell(1), 0, 1);
errorbar_plot(norm_coef_cell(2:7), 1, 2:7);
errorbar_plot(norm_coef_cell(8), 0, 8);
errorbar_plot(norm_coef_cell(9), 0, 9);

% tests
for i = 1:9
    disp(['probe ' num2str(i)])
   [~,coef_norm_p] =  ttest(norm_coef_cell{i});
end
xlim([0.5 num_probes+0.5])
plot(xlim, [1 1].*0, 'k--')
set(gca,'TickLength',[0, 0]); box off;
xlabel('Probe number')
ylabel('Wait times (s)')
title('Normal curve height')
%}



% Intercept of curve
intercept_coef = squeeze(coefficients(:,1,:))';
intercept_coef_cell = cell(1,num_probes);
for iprobe = 1:num_probes
    intercept_coef_cell{iprobe} = intercept_coef(:,iprobe);
end
%{
figure; hold on;
errorbar_plot(intercept_coef_cell(1), 0, 1);
errorbar_plot(intercept_coef_cell(2:7), 1, 2:7);
errorbar_plot(intercept_coef_cell(8), 0, 8);
errorbar_plot(intercept_coef_cell(9), 0, 9);
xlim([0.5 num_probes+0.5])
set(gca,'TickLength',[0, 0]); box off;
xlabel('Probe number')
ylabel('Wait times (s)')
title('Intercept of fit curve')
%}

% Tonal center of curve
center_coef = squeeze(coefficients(:,3,:))';
center_coef_cell = cell(1,num_probes);
for iprobe = 1:num_probes
    center_coef_cell{iprobe} = center_coef(:,iprobe);
end
%{
figure; hold on;
errorbar_plot(center_coef_cell(1), 0, 1);
errorbar_plot(center_coef_cell(2:7), 1, 2:7);
errorbar_plot(center_coef_cell(8), 0, 8);
errorbar_plot(center_coef_cell(9), 0, 9);
xlim([0.5 num_probes+0.5])
set(gca,'TickLength',[0, 0]); box off;
xlabel('Probe number')
ylabel('Curve centers (tone number)')
title('Center of fit curve')
%}



%% Use coefficients to estimate discrimination of every problem tone pair during every probe test

%
% preallocate (probe, problem, subj)
rich_waits_problem = nan(num_probes, num_problems, num_subjects);
poor_waits_problem = nan(num_probes, num_problems, num_subjects);
rich_waits_problem_normal = nan(num_probes, num_problems, num_subjects);
poor_waits_problem_normal = nan(num_probes, num_problems, num_subjects);
rich_waits_problem_logistic = nan(num_probes, num_problems, num_subjects);
poor_waits_problem_logistic = nan(num_probes, num_problems, num_subjects);
diff_waits_problem = nan(num_probes, num_problems, num_subjects);
diff_waits_problem_normal = nan(num_probes, num_problems, num_subjects);
diff_waits_problem_logistic = nan(num_probes, num_problems, num_subjects);

% iterate through all probes
for isubj = 1:size(probe_wait_times,3)
    for iproblem = 1:num_problems
        for iprobe = 1:num_probes
            
            % tone numbers
            rich_tone_num = all_prob_tone_nums(iproblem, 1);
            poor_tone_num = all_prob_tone_nums(iproblem, 2);
            
            
            % tone wait times (during probe) FULL MODEL
            rich_waits_problem(iprobe, iproblem, isubj) = modelFun(coefficients(iprobe,:,isubj), rich_tone_num);
            poor_waits_problem(iprobe, iproblem, isubj) = modelFun(coefficients(iprobe,:,isubj), poor_tone_num);
            
            % tone wait times (during probe) NORMAL DISTRIBUTION ONLY
            normal_only_coefs = coefficients(iprobe,:,isubj);
            normal_only_coefs(3) = realmin;
            
            norm_rich = modelFun(normal_only_coefs, rich_tone_num);
            norm_poor = modelFun(normal_only_coefs, poor_tone_num);
            
            rich_waits_problem_normal(iprobe, iproblem, isubj) = modelFun(normal_only_coefs, rich_tone_num);
            poor_waits_problem_normal(iprobe, iproblem, isubj) = modelFun(normal_only_coefs, poor_tone_num);
            
            % tone wait times (during probe) LOGISTIC DISTRIBUTION ONLY
            logistic_only_coefs = coefficients(iprobe,:,isubj);
            logistic_only_coefs(4) = realmin;
            
            log_rich = modelFun(logistic_only_coefs, rich_tone_num);
            log_poor = modelFun(logistic_only_coefs, poor_tone_num);
            
            rich_waits_problem_logistic(iprobe, iproblem, isubj) = modelFun(logistic_only_coefs, rich_tone_num);
            poor_waits_problem_logistic(iprobe, iproblem, isubj) = modelFun(logistic_only_coefs, poor_tone_num);
            
            
            % difference between rich and poor tone wait times (during probe)
            diff_waits_problem(iprobe, iproblem, isubj)...
                = rich_waits_problem(iprobe, iproblem, isubj) - poor_waits_problem(iprobe, iproblem, isubj);
            diff_waits_problem_normal(iprobe, iproblem, isubj)...
                = rich_waits_problem_normal(iprobe, iproblem, isubj) - poor_waits_problem_normal(iprobe, iproblem, isubj);
            diff_waits_problem_logistic(iprobe, iproblem, isubj)...
                = rich_waits_problem_logistic(iprobe, iproblem, isubj) - poor_waits_problem_logistic(iprobe, iproblem, isubj);

        end
    end
end

% diff_waits_problem_preprobe is prediction score

% plot average difference between each problem tone pair during each probe
figure
ampm_pcolor(nanmean(diff_waits_problem,3))
colorbar
set(gca,'TickLength',[0, 0]); box off;
set(gcf, 'Position', [1092 711 467 628])
title('diff waits problem FULL')
ylabel('Probe number')
xlabel('Problem number')

figure
ampm_pcolor(nanmean(diff_waits_problem_normal,3))
colorbar
set(gca,'TickLength',[0, 0]); box off;
set(gcf, 'Position', [1092 711 467 628])
title('diff waits problem NORMAL')
ylabel('Probe number')
xlabel('Problem number')

figure
ampm_pcolor(nanmean(diff_waits_problem_logistic,3))
colorbar
set(gca,'TickLength',[0, 0]); box off;
set(gcf, 'Position', [1092 711 467 628])
title('diff waits problem LOGISTIC')
ylabel('Probe number')
xlabel('Problem number')


% discrimination of previous problem during each probe
diff_waits_problem_preprobe = nan(num_subjects, num_problems);
diff_waits_problem_postprobe = nan(num_subjects, num_problems);
diff_waits_problem_2postprobe = nan(num_subjects, num_problems-1);
all_other_past_probes = nan(num_subjects, num_problems);
all_future_probes = nan(num_subjects, num_problems);
for isubj = 1:num_subjects
    
    % for current subject
    subj_dwp_preprobe = diff_waits_problem(1:6, :, isubj);
    subj_dwp_postprobe = diff_waits_problem(2:7, :, isubj);
    subj_dwp_2postprobe = diff_waits_problem(3:7, :, isubj);
    
    % extract probe and problem matches
    diff_waits_problem_preprobe(isubj,:) = subj_dwp_preprobe(logical(eye(size(subj_dwp_preprobe))));
    diff_waits_problem_postprobe(isubj,:) = subj_dwp_postprobe(logical(eye(size(subj_dwp_postprobe))));
    diff_waits_problem_2postprobe(isubj,:) = subj_dwp_2postprobe(logical(eye(size(subj_dwp_2postprobe))));
    
    % extract additional probe and problem matches
    for iprobe = 3:7
        all_other_past_probes(isubj, iprobe-1) = nanmean(diff_waits_problem(iprobe, 1:iprobe-2, isubj));
    end
    for iprobe = 1:6
        all_future_probes(isubj, iprobe) = nanmean(diff_waits_problem(iprobe, iprobe:end, isubj));
    end
    
    
end

%
% plot pcolor figure as a line
figure; hold on

    % last problem
    hold_cell = cell(1,size(diff_waits_problem_postprobe,2)); 
    for ic = 1:size(diff_waits_problem_postprobe,2)
        hold_cell{ic} = diff_waits_problem_postprobe(:,ic);
        %[~,last_prob_p] = ttest(diff_waits_problem_postprobe(:,ic))
    end
    errorbar_plot(hold_cell, 1, 2:7, [96 96 96]./255, [0 0 0])

    % all other past problems
    hold_cell = cell(1,size(all_other_past_probes,2)); 
    for ic = 1:size(all_other_past_probes,2)
        hold_cell{ic} = all_other_past_probes(:,ic);
        %[~,past_probs_p] = ttest(all_other_past_probes(:,ic))
    end
    errorbar_plot(hold_cell, 1, 2:7, [238 85 85]./255, [193 24 24]./255)

    % all future problems
    hold_cell = cell(1,size(all_future_probes,2)); 
    for ic = 1:size(all_future_probes,2)
        hold_cell{ic} = all_future_probes(:,ic);
        %[~,future_probs_p] = ttest(all_future_probes(:,ic))
    end
    errorbar_plot(hold_cell, 1, 1:6, [159 188 239]./255, [7 63 159]./255)
    %}

    xlim([0.5 7.5])
    xlabel('Probe number')
    ylabel('Rich - poor tone wait (s)')
    hold on; plot(xlim, [1 1].*0, 'k--')
% individual bars
%{
for ip = 3:6
    figure; 
    errorbar_barplot([{diff_waits_problem_postprobe(:,ip-1)} {all_other_past_probes(:,ip-1)} {all_future_probes(:,ip)}])
    [~,pastVfuture_prob_p] = ttest(all_other_past_probes(:,ip-1), all_future_probes(:,ip))
    title(['Probe number ' num2str(ip) '; past v future pval=' num2str(pastVfuture_prob_p)])
end
%}
        
    
% discrimination of previous problem during each probe (NORMAL COMPONENT)
diff_waits_problem_preprobe_normal = nan(num_subjects, num_problems);
diff_waits_problem_postprobe_normal = nan(num_subjects, num_problems);
for isubj = 1:num_subjects
    
    % for current subject
    subj_dwp_preprobe_normal = diff_waits_problem_normal(1:6, :, isubj);
    subj_dwp_postprobe_normal = diff_waits_problem_normal(2:7, :, isubj);
    
    % extract probe and problem matches
    diff_waits_problem_preprobe_normal(isubj,:) = subj_dwp_preprobe_normal(logical(eye(size(subj_dwp_preprobe_normal))));
    diff_waits_problem_postprobe_normal(isubj,:) = subj_dwp_postprobe_normal(logical(eye(size(subj_dwp_postprobe_normal))));
end

% discrimination of previous problem during each probe (LOGISTIC COMPONENT)
diff_waits_problem_preprobe_logistic = nan(num_subjects, num_problems);
diff_waits_problem_postprobe_logistic = nan(num_subjects, num_problems);
for isubj = 1:num_subjects
    
    % for current subject
    subj_dwp_preprobe_logistic = diff_waits_problem_logistic(1:6, :, isubj);
    subj_dwp_postprobe_logistic = diff_waits_problem_logistic(2:7, :, isubj);
    
    % extract probe and problem matches
    diff_waits_problem_preprobe_logistic(isubj,:) = subj_dwp_preprobe_logistic(logical(eye(size(subj_dwp_preprobe_logistic))));
    diff_waits_problem_postprobe_logistic(isubj,:) = subj_dwp_postprobe_logistic(logical(eye(size(subj_dwp_postprobe_logistic))));
end
%}




%% Coefficient changes
logist_ceof_delta = logist_ceof(:, 2:end) - logist_ceof(:, 1:end-1);
logist_ceof_corrected_delta = logist_ceof_corrected(:, 2:end) - logist_ceof_corrected(:, 1:end-1);
norm_coef_delta = norm_coef(:, 2:end) - norm_coef(:, 1:end-1);
intercept_coef_delta = intercept_coef(:, 2:end) - intercept_coef(:, 1:end-1);
center_coef_delta = center_coef(:, 2:end) - center_coef(:, 1:end-1);




%% probe over probe change (r val)

% number of probes
num_probes = size(probe_wait_times,1);
num_tones = size(probe_wait_times,2);
num_subjects = size(probe_wait_times,3);

% preallocate
pop_rvals = nan(num_subjects, num_probes-1);

% iterate through subjects
for isubj = 1:num_subjects
   for iprobe = 1:num_probes-1
       preprobe_waits = probe_wait_times(iprobe,:,isubj)';
       postprobe_waits = probe_wait_times(iprobe+1,:,isubj)';
       nnan_idx = ~isnan(preprobe_waits) & ~isnan(postprobe_waits);
       pop_rvals(isubj, iprobe) = corr(preprobe_waits(nnan_idx), postprobe_waits(nnan_idx));
   end
end
figure; ampm_pcolor(pop_rvals); title('pop rvals')



%% Population remapping


% neuronal population overlap
overlap_mtx_3d = overlapping_population_average(unq_subjs);

% probe days only
probe_days = [1:3:19 20];
a_idx = zeros(8,8);
a_idx([[zeros(7,1) eye(7,7)];zeros(1,8)]==1) = 1;

overlap_mtx = nan(7, size(overlap_mtx_3d,3));
for isubj = 1:size(overlap_mtx_3d,3)
    overlap_mtx_3d_subj = overlap_mtx_3d(:,:,isubj);
    overlap_mtx_3d_subj = overlap_mtx_3d_subj(probe_days,probe_days);
    overlap_mtx(:,isubj) = overlap_mtx_3d_subj(logical(a_idx)); 
end


% mean correlation between event-aligned firing-rate tuning curves
% see cellCorrs_v_probeCorrs
load('cellCorrs_v_probeCorrs_prep.mat', 'subj_corr_cell_mevar', 'subj_corr_cell_hivar')
if contains(training_group, 'mevar')
    cell_corrs = subj_corr_cell_mevar;
elseif contains(training_group, 'hivar')
    cell_corrs = subj_corr_cell_hivar;
end
            
% minimum number of cells
min_cell = 3;
% means only
cell_corr_means = nan(length(cell_corrs{1}),length(cell_corrs{1}),length(cell_corrs));
for isubj = 1:length(cell_corrs)
    for isesh1 = 1:size(cell_corrs{isubj},1)
        for isesh2 = 1:size(cell_corrs{isubj},2)
            
            % no redundancy
            if isesh2>=isesh1
                continue
            end
            if length(cell_corrs{isubj}{isesh1, isesh2})>=min_cell
                %cell_corr_means(isesh1,isesh2,isubj) = nanmean(atanh(cell_corrs{isubj}{isesh1, isesh2}));
                cell_corr_means(isesh1,isesh2,isubj) = nanmean((cell_corrs{isubj}{isesh1, isesh2}));
            else
                cell_corr_means(isesh1,isesh2,isubj) = nan;
            end
        end
    end
end


% test plot
cell_corr_means
for isubj = 1:size(cell_corr_means,3)
figure; imagesc(cell_corr_means(:,:,isubj))
end

% probe days only
probe_days = [1:3:19 20];
cell_corr_means = cell_corr_means(probe_days,probe_days,:)

% just pre-to-post probe remaps
a_idx = zeros(8,8);
a_idx([[zeros(7,1) eye(7,7)];zeros(1,8)]==1) = 1;
cell_corr_means_mtx = nan(7,size(cell_corr_means,3));
for i = 1:size(cell_corr_means,3)
    cell_corr_means_hold = cell_corr_means(:,:,i)'; % transpose to deal with the no redundancy decision above
    cell_corr_means_mtx(:,i) = cell_corr_means_hold(logical(a_idx)); 
end
cell_corr_means_mtx



%% Questions

    
    % quality of prediction related to REMAPPING
    
    % population overlap
    accuracy_of_prediction = diff_waits_problem_preprobe;
    pop_overlap = overlap_mtx(1:6,:)';
        
    title_string = [remove_underscore(training_group) ': Preprobe problem discrim predicts probe-to-probe population overlap'];  
    q7 = lme_function(model_str_against, pop_overlap, accuracy_of_prediction, subj_num_mtx, stage_num_mtx);
    plot2var_corr(accuracy_of_prediction, pop_overlap, title_string)
    plot2var(accuracy_of_prediction, pop_overlap, title_string)    
    
    
    
    % remapping
    accuracy_of_prediction = diff_waits_problem_preprobe;
    remap_score = cell_corr_means_mtx(1:6,:)';
    
        
    
    size(accuracy_of_prediction)
    size(remap_score)
    
    title_string = [remove_underscore(training_group) ': Preprobe problem discrim predicts probe-to-probe population firing code stability'];
    q8 = lme_function(model_str_against, remap_score, accuracy_of_prediction, subj_num_mtx, stage_num_mtx);
    plot2var_corr(accuracy_of_prediction, remap_score, title_string)
    plot2var(accuracy_of_prediction, remap_score, title_string) 

    
    
    
    
    
 

end


% Figure functions
function plot2var(var1, var2, title_str)
% vars need to be matrices with (subjects, samples)

    % colormap
    colors = bone(10); colors = colors(2:9,:);

    % open figure
    figure;
    
    % mean centered correlation dot plot
    subplot(1,2,1); hold on
    for iprobe = 1:size(var1,2)

        % mean center
        %
        var1_mc(:, iprobe) = var1(:, iprobe) - nanmean(var1(:, iprobe));
        var2_mc(:, iprobe) = var2(:, iprobe) - nanmean(var2(:, iprobe));
        %}
        plot(var1_mc(:, iprobe), var2_mc(:, iprobe), 'o', 'color', colors(iprobe,:))
    end
    [r, p] = fit_line(var1_mc(:), var2_mc(:), 0);
    set(gca,'TickLength',[0, 0]); box off;
    title(['r=' num2str(r) ', p=' num2str(p)])
    xlabel([remove_underscore(inputname(1)) ', mean centered'])
    ylabel([remove_underscore(inputname(2)) ', mean centered'])
    axis square
    hold on; plot(xlim, [0 0], 'k--')
    hold on; plot([0 0], ylim, 'k--')


    % dot plot with lines
    subplot(1,2,2); hold on
    for isubj = 1:size(var1,1)
        for iprobe = 1:size(var1,2)
            if iprobe < size(var1,2)
                plot(var1(isubj, [iprobe iprobe+1]), var2(isubj, [iprobe iprobe+1]), '-', 'color', colors(iprobe,:), 'linewidth', 0.5)
            end
            plot(var1(isubj, iprobe), var2(isubj, iprobe), 'o', 'color', colors(iprobe,:))
        end
    end
    for iprobe = 1:size(var1,2)
        if iprobe < size(var1,2)
            plot(nanmean(var1(:, [iprobe iprobe+1])), nanmean(var2(:, [iprobe iprobe+1])), '-', 'color', colors(iprobe,:), 'linewidth', 6)
        end
        plot(nanmean(var1(:, iprobe)), nanmean(var2(:, iprobe)), '.', 'color', colors(iprobe,:), 'markersize', 60)
    end
    [r, p] = fit_line(var1(:), var2(:), 0);
    title(['r=' num2str(r) ', p=' num2str(p)])
    set(gca,'TickLength',[0, 0]); box off;
    xlabel(remove_underscore(inputname(1)))
    ylabel(remove_underscore(inputname(2)))
    axis square
    hold on; plot(xlim, [0 0], 'k--')
    hold on; plot([0 0], ylim, 'k--')
    
    % overarching title
    sgtitle(title_str)
end

% plot correlations for each probe/problem
function plot2var_corr(var1, var2, title_str)

    % if nan in one, nan in both
    nan_idx = isnan(var1) | isnan(var2);
    var1(nan_idx) = nan;
    var2(nan_idx) = nan;



    figure
    for iprob = 1:size(var1,2)
        
        if all(isnan(var1(:,iprob))) || all(isnan(var2(:,iprob)))
             continue
        end
        
        subplot(1,size(var1,2),iprob) 
        hold on
        
        [r,p] = fit_line(var1(:,iprob), var2(:,iprob));
        axis square
        plot(xlim, [1 1].*0, 'k--')
        plot([1 1].*0, ylim, 'k--')
        xlabel(remove_underscore(inputname(1)))
        ylabel(remove_underscore(inputname(2)))
        title(['Prob ' num2str(iprob) '; r=' num2str(r) ', p=' num2str(p)])
    end
    sgtitle(title_str)
    set(gcf, 'Position', [680 558 1533 420])
end


function lme = lme_function(model_str, var1, var2, subj_num_mtx, stage_num_mtx)
    % function testing whether var2 (input 3) predicts var1 (input 2)
    
    
    % if nan in one, nan in both
 %   nan_idx = isnan(var1) | isnan(var2);
 %   var1(nan_idx) = nan_idx;
 %   var2(nan_idx) = nan_idx;

    
    
    % update model string with input names

        %var1
        var1_str_starts = strfind(model_str, 'var1');
        for v1ss = fliplr(var1_str_starts)
            model_str = [model_str(1:v1ss-1) inputname(2) model_str((v1ss+length('var1')):end)];
        end

        %var2
        var2_str_starts = strfind(model_str, 'var2');
        for v2ss = fliplr(var2_str_starts)
            model_str = [model_str(1:v2ss-1) inputname(3) model_str((v2ss+length('var2')):end)];
        end
        
    % mixed model
    tbl = table(var1(:), var2(:), subj_num_mtx(:), stage_num_mtx(:),...
        'VariableNames',{inputname(2), inputname(3), 'subject', 'problem'});
    tbl.subject = categorical(tbl.subject);
    %tbl.problem = categorical(tbl.problem);
    lme = fitlme(tbl, model_str);
                
end

function lme = lme_function_double(model_str, var1, var2, var3, subj_num_mtx, stage_num_mtx)
    % function testing whether var2 (input 3) predicts var1 (input 2)
    
    % update model string with input names

        %var1
        var1_str_starts = strfind(model_str, 'var1');
        for v1ss = fliplr(var1_str_starts)
            model_str = [model_str(1:v1ss-1) inputname(2) model_str((v1ss+length('var1')):end)];
        end

        %var2
        var2_str_starts = strfind(model_str, 'var2');
        for v2ss = fliplr(var2_str_starts)
            model_str = [model_str(1:v2ss-1) inputname(3) model_str((v2ss+length('var2')):end)];
        end
        
        %var3
        var3_str_starts = strfind(model_str, 'var3');
        for v3ss = fliplr(var3_str_starts)
            model_str = [model_str(1:v3ss-1) inputname(4) model_str((v3ss+length('var3')):end)];
        end
        
    % mixed model
    tbl = table(var1(:), var2(:), var3(:), subj_num_mtx(:), stage_num_mtx(:),...
        'VariableNames',{inputname(2), inputname(3), inputname(4), 'subject', 'problem'});
    tbl.subject = categorical(tbl.subject);
    tbl.problem = categorical(tbl.problem);
    lme = fitlme(tbl, model_str);
                
end





function [logist_coef_cell_corrected, norm_coef_cell, diff_waits_problem_preprobe, diff_waits_problem_postprobe, correlation_plot_cell, diff_waits_problem] = ALL_probe_vs_behavior_fit(training_group)
% Used to compare training problem behavior with probe behavior in
% every which way. produces many plots. additional ones found near the
% bottom

% data of interest out
% correlation_plot_cell{1} = first day discrim
% correlation_plot_cell{2} = days to crit
% correlation_plot_cell{3} = amount learning
% correlation_plot_cell{4} = probe-over-probe similarity

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
for iplot = 1:7; close; end % close legacy plots
% firstLast_waitTimes_cell:
%   6 cells (1 per problem)
%       2 cells per subj (rows: subj, columns: firstDay, lastDay)
%           2 cells (rich, poor)
%               vector of wait times

% first and last day dprimes (subj, problem)
first_day_d = squeeze(firstLast_dprime_mtx(:,1,:));
last_day_d = squeeze(firstLast_dprime_mtx(:,2,:));

% learning over first two days
%{
TwoDay_learn_mtx = nan(size(all_TwoDay_learn_mtx,1), num_problems);
for iprob = 1:num_problems
    TwoDay_learn_mtx(:,iprob) = all_TwoDay_learn_mtx(:, 2, iprob)-all_TwoDay_learn_mtx(:, 1, iprob);
end
%}

%% all probe wait times

% wait times from every probe for every subject (probe,tone,subj)
[probe_wait_times] = ALL_probe_wait_times(training_group, 1); % model fit, interp & extrap, no smooth

% number of probes
num_probes = size(probe_wait_times,1);
num_tones = size(probe_wait_times,2);
num_subjects = size(probe_wait_times,3);

% zscored probe wait times
%{
probe_wait_times_z = nan(num_probes, num_tones, num_subjects);
for iprobe = 1:num_probes
    for isubj = 1:num_subjects
        probe_wait_times_z(iprobe, :, isubj) = zscore_mtx(probe_wait_times(iprobe, :, isubj)')';
    end
end
%}



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

        % compute coefficients
        try
            mean_wait_times = nanmean(probe_wait_times(iprobe,:,isubj));
            %[~, coefEsts, modelFun] = ampm_normal_logistic_fit_algo(1:num_tones, wait_times, [mean_wait_times 20 0 1]);
            [~, coefEsts, modelFun] = ampm_normal_logistic_fit(1:num_tones, wait_times, [mean_wait_times 20 0 1]);
  
        % if modelling fails, set to mean
        catch 
            try
                mean_wait_times = nanmean(probe_wait_times(iprobe,:,isubj));
                [~, coefEsts, modelFun] = ampm_normal_logistic_fit_algo(1:num_tones, wait_times, [mean_wait_times-5 21 0 10]);
            catch
                
                try
                    mean_wait_times = nanmean(probe_wait_times(iprobe,:,isubj));
                    [~, coefEsts, modelFun] = ampm_normal_logistic_fit(1:num_tones, wait_times, [mean_wait_times 20 0 1]);
                
                catch
                    figure; plot(1:num_tones, wait_times,'o'); title(['subj ' num2str(isubj) '; probe ' num2str(iprobe)])
                    coefEsts = [nanmean(probe_wait_times(iprobe,:,isubj)) zeros(1,num_coefs-1)];
                end
            end
        end
        
        % load coefficients
        coefficients(iprobe,:,isubj) = coefEsts;
        
        % santity check coefs
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
    end
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
%
figure; hold on;
errorbar_plot(logist_coef_cell_corrected(1), 0, 1);
errorbar_plot(logist_coef_cell_corrected(2:7), 1, 2:7);
errorbar_plot(logist_coef_cell_corrected(8), 0, 8);
if num_probes == 9
    errorbar_plot(logist_coef_cell_corrected(9), 0, 9);
end

%tests
for i = 1:num_probes
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
%
figure; hold on;
errorbar_plot(norm_coef_cell(1), 0, 1);
errorbar_plot(norm_coef_cell(2:7), 1, 2:7);
errorbar_plot(norm_coef_cell(8), 0, 8);
if num_probes ==9
errorbar_plot(norm_coef_cell(9), 0, 9);
end

% tests
%{
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


% Tonal center of curve
center_coef = squeeze(coefficients(:,3,:))';
center_coef_cell = cell(1,num_probes);
for iprobe = 1:num_probes
    center_coef_cell{iprobe} = center_coef(:,iprobe);
end




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
            
            rich_waits_problem_normal(iprobe, iproblem, isubj) = modelFun(normal_only_coefs, rich_tone_num);
            poor_waits_problem_normal(iprobe, iproblem, isubj) = modelFun(normal_only_coefs, poor_tone_num);
            
            % tone wait times (during probe) LOGISTIC DISTRIBUTION ONLY
            logistic_only_coefs = coefficients(iprobe,:,isubj);
            logistic_only_coefs(4) = realmin;
            
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

    % previous problem
    %
    figure; hold on
    hold_cell = cell(1,size(diff_waits_problem_postprobe,2)); 
    for ic = 1:size(diff_waits_problem_postprobe,2)
        hold_cell{ic} = diff_waits_problem_postprobe(:,ic);
    end
    errorbar_plot(hold_cell, 1, 2:7, [96 96 96]./255, [0 0 0])
    xlim([0.5 7.5]); ylim([-30 30])
    xlabel('Probe number')
    ylabel('Rich - poor tone wait (s)')
    hold on; plot(xlim, [1 1].*0, 'k--')
    title('Previous problem')
    %}

    % immediate future problem
    %
    figure; hold on
    hold_cell = cell(1,size(diff_waits_problem_preprobe,2)); 
    for ic = 1:size(diff_waits_problem_preprobe,2)
        hold_cell{ic} = diff_waits_problem_preprobe(:,ic);
    end
    errorbar_plot(hold_cell, 1, 1:6, [0 128 0]./255, [34 139 34]./255)
    %}


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



%% Questions

    
% quality of prediction related to first day performance

    accuracy_of_prediction = diff_waits_problem_preprobe;
    first_day_discrimination = first_day_d;
    correlation_plot_cell{1} = first_day_discrimination;
    title_string = [remove_underscore(training_group) ': Preprobe problem discrim predicts first day performance'];
    q1 = lme_function(model_str_against, first_day_discrimination, accuracy_of_prediction, subj_num_mtx, stage_num_mtx);
    plot2var(accuracy_of_prediction, first_day_discrimination, title_string)
    
    
% quality of prediction related to days to criterion

    accuracy_of_prediction = diff_waits_problem_preprobe;
    days_to_criterion = days_to_crit;
    correlation_plot_cell{2} = days_to_criterion;
    title_string = [remove_underscore(training_group) ': Preprobe problem discrim predicts first day performance'];
    q2 = lme_function(model_str_against, days_to_criterion, accuracy_of_prediction, subj_num_mtx, stage_num_mtx);
    plot2var(accuracy_of_prediction, days_to_criterion, title_string)
 
    
% quality of prediction related to post-probe problem discrim
    
    accuracy_of_prediction = diff_waits_problem_preprobe;
    amount_learning = last_day_d - first_day_d;
    correlation_plot_cell{3} = amount_learning;
    title_string = [remove_underscore(training_group) ': Preprobe problem discrim predicts amount learning'];
    q3 = lme_function(model_str_against, amount_learning, accuracy_of_prediction, subj_num_mtx, stage_num_mtx);
    plot2var(accuracy_of_prediction, amount_learning, title_string)    

    
% quality of prediction related to probe-over-probe change

    accuracy_of_prediction = diff_waits_problem_preprobe;
    probe_over_probe_corr = pop_rvals(:,1:6);
    correlation_plot_cell{4} = probe_over_probe_corr;
    title_string = [remove_underscore(training_group) ': Preprobe problem discrim predicts probe-over-probe similarity'];
    q4 = lme_function(model_str_against, probe_over_probe_corr, accuracy_of_prediction, subj_num_mtx, stage_num_mtx);
    plot2var(accuracy_of_prediction, probe_over_probe_corr, title_string)
    

% quality of prediction related to last day performance

    accuracy_of_prediction = diff_waits_problem_preprobe;
    last_day_discrimination = last_day_d;
    correlation_plot_cell{5} = last_day_discrimination;
    title_string = [remove_underscore(training_group) ': Preprobe problem discrim predicts last day performance'];
    q1 = lme_function(model_str_against, last_day_discrimination, accuracy_of_prediction, subj_num_mtx, stage_num_mtx);
    plot2var(accuracy_of_prediction, last_day_discrimination, title_string)
    
    

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

function lme = lme_function(model_str, var1, var2, subj_num_mtx, stage_num_mtx)
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
        
    % mixed model
    tbl = table(var1(:), var2(:), subj_num_mtx(:), stage_num_mtx(:),...
        'VariableNames',{inputname(2), inputname(3), 'subject', 'problem'});
    tbl.subject = categorical(tbl.subject);
    lme = fitlme(tbl, model_str);
                
end






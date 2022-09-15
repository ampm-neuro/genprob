function ALL_probe_vs_behavior_new(training_group)
% Used to compare training problem behavior with probe behavior in
% every which way. produces many plots. additional ones found near the
% bottom



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
        all_prob_tone_nums(iprob, itone) = find(ismember(unqfrq41, all_prob_tones(iprob,itone)));
    end
end



%% Discrimination problem behavior

% all wait times over learning
[firstLast_dprime_mtx, firstLast_waitTimes_cell, days_to_crit, all_TwoDay_learn_mtx] = ALL_stage_learn(['two_tone\' training_group], 1:6, 2);
% firstLast_waitTimes_cell:
%   6 cells (1 per problem)
%       2 cells per subj (rows: subj, columns: firstDay, lastDay)
%           2 cells (rich, poor)
%               vector of wait times

% first and last day dprimes (subj, problem)
first_day_d = squeeze(firstLast_dprime_mtx(:,1,:));
last_day_d = squeeze(firstLast_dprime_mtx(:,2,:));

% learning over first two days
TwoDay_learn_mtx = nan(size(all_TwoDay_learn_mtx,1), num_problems);
for iprob = 1:num_problems
    TwoDay_learn_mtx(:,iprob) = all_TwoDay_learn_mtx(:, 2, iprob)-all_TwoDay_learn_mtx(:, 1, iprob);
end



%% all probe wait times

% smooth wait times from every probe for every subject (probe,tone,subj)
[probe_wait_times, unq_subjects_probe_all] = ALL_probe_wait_times(training_group, 2);

% number of probes
num_probes = size(probe_wait_times,1);
num_tones = size(probe_wait_times,2);
num_subjects = size(probe_wait_times,3);

% zscored probe wait times
probe_wait_times_z = nan(num_probes, num_tones, num_subjects);
for iprobe = 1:num_probes
    for isubj = 1:num_subjects
        probe_wait_times_z(iprobe,:,isubj) = zscore_mtx(probe_wait_times(iprobe,:,isubj)')';
    end
end



%% Colors

% colors
colors = bone(10); colors = colors(2:num_probes,:);



%% Probe wait time plots

% iterate through probes
for iprobe = 1:num_probes

    % wait times on this probe from all subjects (subj, tone)
    probes_local = squeeze(probe_wait_times(iprobe, :, :))';
    
    % plot average wait times at each frequency
    %
    figure; hold on
    for isubj = 1:num_subjects 
        %plot(probes_local(isubj,:),'color', [102 178 204]./255)
    end
    errorbar(nanmean(probes_local,1), nanstd(probes_local,[],1)./sqrt(sum(~isnan(probes_local),1)), 'linewidth', 2, 'color', [0 124 204]./255);
    %ylim([0 40])
    ylim([0 25])
    ylabel('Wait time (s)')
    xlim([-.5 42.5])
    xlabel('Tone')
    set(gca,'TickLength',[0, 0]); box off;
    title(['Probe ' num2str(iprobe)])
    %{
    if iprobe >=2 && iprobe <=7
        hold on
        rich_bounds_prob('hivar', iprobe-1, 1)
    end
    var_name = ['probe0' num2str(iprobe)]; print(['C:\Users\ampm1\Documents\LabMeet\LabMeet07\hivar_probes\' var_name], '-dpdf', '-painters', '-bestfit')
    %}
end


%% probe to probe similarity

% preallocate probe behavior
probe_delta_r = nan(num_subjects, num_probes-1);

% iterate through probes
for iprobe = 1:num_probes-1

    % wait times on this probe from all subjects (subj, tone)
    probes_pre = squeeze(probe_wait_times(iprobe, :, :))';
    probes_post = squeeze(probe_wait_times(iprobe+1, :, :))';

    % compute similarity measures(correlation and difference) for all subjects
    for isubjCorr = 1:num_subjects
        probe_delta_r(isubjCorr, iprobe) = corr(probes_pre(isubjCorr, :)', probes_post(isubjCorr, :)');
    end
end


%% Plot mean probe over probe change in terms of similarlity (correlation)

% gather and transform r values
pd_full = atanh(probe_delta_r); %fischer z transform for normality

% prepare r values for errorbar_plot
pd_full_vect = [];
for iprobe = 1:num_probes-2
    pd_full_vect = [pd_full_vect {pd_full(:,iprobe)}];
    
    [a b c d] = ttest(pd_full(:,1), pd_full(:,iprobe))
    
end

    % plot r values across learning in an lineplot with errorbars
    figure; hold on
    errorbar_plot(pd_full_vect(1:6), 1, 1:6)
    errorbar_plot(pd_full_vect(7), 0, 7)
    xlim([.1 7.9])
    title('mean probe over probe change')
    xlabel('Training problem number')
    ylabel('Probe over probe similarity (r)')
    hold on; plot(xlim, [0 0], 'k--')
    ylim([-1 1])
    maxr = .975; tic_vect = [-maxr -.9 -.8 -.6 -.3 0 .3 .6 .8 .9 maxr];
    ylim(atanh([-maxr maxr])); yticks(atanh(tic_vect)); yticklabels(tic_vect)

    % plot just first and last r values in a barplot
    figure
    errorbar_barplot([pd_full_vect(1) pd_full_vect(end)], 1);
    [~, r_FirstLast] = ttest(pd_full_vect{1}, pd_full_vect{end});
    title(['Similarity: first two and last two probes ; p = ' num2str(r_FirstLast)])
    ylabel('Probe over probe similarity (r)')
    ylim([-1 1])
    ylim(atanh([-maxr maxr])); yticks(atanh(tic_vect)); yticklabels(tic_vect)



%% Waits (z) to each problem tone DURING PROBES

% preallocate (probe, problem, subj)
rich_waits_problem = nan(num_probes, num_problems, num_subjects);
poor_waits_problem = nan(num_probes, num_problems, num_subjects);
diff_waits_problem = nan(num_probes, num_problems, num_subjects);
rich_zwaits_problem = nan(num_probes, num_problems, num_subjects);
poor_zwaits_problem = nan(num_probes, num_problems, num_subjects);
zdiff_waits_problem = nan(num_probes, num_problems, num_subjects);

% iterate through all probes
for isubj = 1:size(probe_wait_times,3)
    for iproblem = 1:num_problems
        for iprobe = 1:num_probes
            
            % tone numbers
            rich_tone_num = all_prob_tone_nums(iproblem, 1);
            poor_tone_num = all_prob_tone_nums(iproblem, 2);
            
            
            % tone wait times (during probe)
            rich_waits_problem(iprobe, iproblem, isubj) = probe_wait_times(iprobe, rich_tone_num, isubj);
            poor_waits_problem(iprobe, iproblem, isubj) = probe_wait_times(iprobe, poor_tone_num, isubj);
            
            % difference between rich and poor tone wait times (during probe)
            diff_waits_problem(iprobe, iproblem, isubj)...
                = rich_waits_problem(iprobe, iproblem, isubj) - poor_waits_problem(iprobe, iproblem, isubj);

            
            % tone zscores (during probe)
            rich_zwaits_problem(iprobe, iproblem, isubj) = probe_wait_times_z(iprobe, rich_tone_num, isubj);
            poor_zwaits_problem(iprobe, iproblem, isubj) = probe_wait_times_z(iprobe, poor_tone_num, isubj);
            
            % zscore difference between rich and poor tone (during probe)
            zdiff_waits_problem(iprobe, iproblem, isubj)...
                = rich_zwaits_problem(iprobe, iproblem, isubj) - poor_zwaits_problem(iprobe, iproblem, isubj);
    
        end
    end
end

% discrimination of previous problem during each probe
diff_waits_mean_problem_preprobe = nan(num_subjects, num_problems);
diff_waits_mean_problem_postprobe = nan(num_subjects, num_problems);
for isubj = 1:num_subjects
    
    % for current subject
    subj_dwp_preprobe = diff_waits_problem(1:6, :, isubj);
    subj_dwp_postprobe = diff_waits_problem(2:7, :, isubj);
    
    % extract probe and problem matches
    diff_waits_mean_problem_preprobe(isubj,:) = subj_dwp_preprobe(logical(eye(size(subj_dwp_preprobe))));
    diff_waits_mean_problem_postprobe(isubj,:) = subj_dwp_postprobe(logical(eye(size(subj_dwp_postprobe))));
end


% average discrimination of every problem during each probe
diff_waits_mean_richVpoor_cell = cell(1,num_probes);
for iprobe = 1:num_probes
    diff_waits_mean_richVpoor_cell{iprobe} = squeeze(mean(diff_waits_problem(iprobe, :, :),2));
end
diff_waits_mean_richVpoor = cell2mat(diff_waits_mean_richVpoor_cell);
diff_waits_mean_richVpoor_preprobe = diff_waits_mean_richVpoor(:,1:6);
diff_waits_mean_richVpoor_postprobe = diff_waits_mean_richVpoor(:,2:7);


% plot wait times at problem tones during probes
%
% line errorbar plot
    all_rich_changes = cell(1, num_probes);
    all_poor_changes = cell(1, num_probes);
    for iprobe = 1:num_probes
        all_rich_changes{iprobe} = squeeze(mean(rich_waits_problem(iprobe, :, :),2)); % average wait to all rich tones (incl future)
        all_poor_changes{iprobe} = squeeze(mean(poor_waits_problem(iprobe, :, :),2)); % average discrim to all poor tones (incl future)
    end
    figure; hold on;
    title('Rich tone and poor tone wait times')

    errorbar_plot(all_rich_changes(1), 1, 1, [0 0.447 0.741]);
    errorbar_plot(all_rich_changes(2:7), 1, 2:7, [0 0.447 0.741]);
    errorbar_plot(all_rich_changes(8), 1, 8, [0 0.447 0.741]);

    errorbar_plot(all_poor_changes(1), 1, 1, [0.85 0.325 0.098]);
    errorbar_plot(all_poor_changes(2:7), 1, 2:7, [0.85 0.325 0.098]);
    errorbar_plot(all_poor_changes(8), 1, 8, [0.85 0.325 0.098]);
    
    xlim([0.5 8.5])
    hold on; plot(xlim, [1 1].*0, 'k--')


    
% plot average discrimination of every problem during each probe
%
    % plot average discrim index
    figure; hold on
    errorbar_plot(diff_waits_mean_richVpoor_cell(1), 1, 1);
    errorbar_plot(diff_waits_mean_richVpoor_cell(2:7), 1, 2:7);
    errorbar_plot(diff_waits_mean_richVpoor_cell(8), 1, 8);
    xlim([0.5 8.5])
    plot(xlim, [1 1].*0, 'k--')
    xlabel('Probe number')
    ylabel('Rich waits - poor waits')
    title('Difference between rich and poor tone responses during probe sessions')
    



%% Waits (z) OUTSIDE each problem tones DURING PROBES

% preallocate (probe, subj)
high_waits = nan(num_probes, num_subjects);
low_waits = nan(num_probes, num_subjects);
diff_waits_outside = nan(num_probes, num_subjects);
high_zwaits = nan(num_probes, num_subjects);
low_zwaits = nan(num_probes, num_subjects);
zdiff_waits_outside = nan(num_probes, num_subjects);

% iterate through all probes
target_tone = [1 num_tones];
for isubj = 1:size(probe_wait_times,3)
    for iprobe = 1:num_probes

        
        % flip target tones according to relative position of rich/poor tones
        %{
        if iprobe > 1 && iprobe < 8
            if all_prob_tone_nums(iprobe-1, 1) > all_prob_tone_nums(iprobe-1, 2)
                target_tone = [num_tones 1];
            else
                target_tone = [1 num_tones];
            end
        else
            target_tone = [1 num_tones];
        end
        %}
        
        % first and last tone waits (high and low may also be richside and
        % poorside depending on above code)
        low_waits(iprobe, isubj) = probe_wait_times(iprobe, target_tone(1), isubj);
        high_waits(iprobe, isubj) = probe_wait_times(iprobe, target_tone(2), isubj);
            
        % zscore difference between outsides of rich and poor tones (during probe)
        diff_waits_outside(iprobe, isubj)...
            = low_waits(iprobe, isubj) - high_waits(iprobe, isubj);
        
        % first and last tone z waits
        low_zwaits(iprobe, isubj) = probe_wait_times_z(iprobe, target_tone(1), isubj);
        high_zwaits(iprobe, isubj) = probe_wait_times_z(iprobe, target_tone(2), isubj);
        
        % zscore difference between outsides of rich and poor tones (during probe)
        zdiff_waits_outside(iprobe, isubj)...
            = low_zwaits(iprobe, isubj) - high_zwaits(iprobe, isubj);
        
    end
end


% plot wait outside problem tones    
%
% line errorbar plot
    all_high_waits = cell(1, num_probes);
    all_low_waits = cell(1, num_probes);
    for iprobe = 1:num_probes
        all_high_waits{iprobe} = high_zwaits(iprobe, :)';
        all_low_waits{iprobe} = low_zwaits(iprobe, :)';
    end
    figure; hold on;
    
    errorbar_plot(all_high_waits(1), 1, 1, [0 0.447 0.741]);
    errorbar_plot(all_high_waits(2:7), 1, 2:7, [0 0.447 0.741]);
    errorbar_plot(all_high_waits(8), 1, 8, [0 0.447 0.741]);

    errorbar_plot(all_low_waits(1), 1, 1, [0.85 0.325 0.098]);
    errorbar_plot(all_low_waits(2:7), 1, 2:7, [0.85 0.325 0.098]);
    errorbar_plot(all_low_waits(8), 1, 8, [0.85 0.325 0.098]);
    
    xlim([0.5 8.5])
    plot(xlim, [1 1].*0, 'k--')
    title('High and low tone responses')


% plot average discrimination of outsides of each problem during each probe
diff_waits_mean_outside = cell(1,num_probes);
for iprobe = 1:num_probes
    diff_waits_mean_outside{iprobe} = zdiff_waits_outside(iprobe,:)';
end

    % plot average discrim index
    figure; hold on
    errorbar_plot(diff_waits_mean_outside(1), 1, 1);
    errorbar_plot(diff_waits_mean_outside(2:7), 1, 2:7);
    errorbar_plot(diff_waits_mean_outside(8), 1, 8);
    xlim([0.5 8.5])
    plot(xlim, [1 1].*0, 'k--')
    xlabel('Probe number')
    ylabel('Low freq end - high freq end')
    title('Difference between high and low tone responses during probe sessions')

    
% CHANGE BEYOND / OUTSIDE OF TONES    
% does the change in dprime from first to last day predict the rich/poor 
% probe side difference

    % discrimination of outside tones on next problem (problems 1:6)  
    % during current probe (probes 1:6)
    probe_discr_next_prob_outside = zdiff_waits_outside(1:6,:)';
    
    % discrimination of outside tones on prior problem (problems 1:6)  
    % during current probe (probes 2:7)
    probe_discr_past_prob_outside = zdiff_waits_outside(2:7,:)';
    
    % change in discrimination of each problem from one probe to the next
    probe_outsidetoneDiscrim_delta = probe_discr_past_prob_outside - probe_discr_next_prob_outside;    

   
    
%% mixed model prep
subj_num_mtx = repmat((1:num_subjects)', 1, num_problems); % subject matrix
stage_num_mtx = repmat(1:6, num_subjects, 1); % stage matrix
model_str = 'var1~var2+(1|subject)+(1|problem)';
model_str_against = 'var1~var2+problem+(1|subject)';


 %% Questions

% does dissagreement between prior probe and new problem predict magnitude
% of probe-over-probe change



% does total learning during problem predict the predict magnitude of 
% probe-over-probe change USE UPCOMING PROBLEM DISCRIM, NOT RichVPoor

    % logistic
    %{
    total_learning = last_day_d; %
    post_probe_outside = probe_discr_past_prob_outside; % post probe only
    title_string = [remove_underscore(training_group) ': Last day predict outside discrim of next probe'];
    plot2var_corr(total_learning, post_probe_outside, title_string)
    plot2var(total_learning, post_probe_outside, title_string)
    %}
    %{
    total_learning = last_day_d-first_day_d; %
    post_probe_outside = probe_discr_past_prob_outside; % post probe only
    title_string = [remove_underscore(training_group) ': Total learning predict outside discrim of next probe'];
    plot2var_corr(total_learning, post_probe_outside, title_string)
    plot2var(total_learning, post_probe_outside, title_string)
    %}
    
    %{
    total_learning = last_day_d-first_day_d; %
    probe_over_probe_outside = abs(probe_outsidetoneDiscrim_delta(:, 1:6)); % delta
    title_string = [remove_underscore(training_group) ': Total learning predict outside discrim delta'];
    %plot2var_corr(total_learning, probe_over_probe_outside, title_string)
    plot2var(total_learning, probe_over_probe_outside, title_string)
    %}
    
    
    % normal
    %{
    total_learning = last_day_d-first_day_d; %
    post_probe_problemDiscrim = diff_waits_mean_problem_postprobe; % post probe only
    title_string = [remove_underscore(training_group) ': Total learning predict probe problem discrim during next probe'];
    plot2var_corr(total_learning, post_probe_problemDiscrim, title_string)
    plot2var(total_learning, post_probe_problemDiscrim, title_string)
    %}
    total_learning = last_day_d-first_day_d; %
    probe_over_probe_problemDiscrim_delta = diff_waits_mean_problem_postprobe - diff_waits_mean_problem_preprobe; % delta
    title_string = [remove_underscore(training_group) ': Total learning predict probe problem discrim delta'];
    %plot2var_corr(total_learning, probe_over_probe_problemDiscrim_delta, title_string)
    plot2var(total_learning, probe_over_probe_problemDiscrim_delta, title_string)
    total_learn_full_pred_probDiscrimDelta = lme_function(model_str, probe_over_probe_problemDiscrim_delta, total_learning, subj_num_mtx, stage_num_mtx);
    total_learn_full_pred_probDiscrimDelta_against = lme_function(model_str_against, probe_over_probe_problemDiscrim_delta, total_learning, subj_num_mtx, stage_num_mtx)
    

% does last day problem discrimination predict post-probe problem discrimination?
%
    
    % normal
    %
    postprobe_problemDiscrim = diff_waits_mean_problem_postprobe; % pre probe only
    last_day_discrimination = last_day_d;
    title_string = [remove_underscore(training_group) ': Last day discrimination predicts post probe discrimination '];
    %plot2var_corr(preprobe_problemDiscrim, first_day_discrimination, title_string)
    plot2var(last_day_discrimination, postprobe_problemDiscrim, title_string)
    last_day_discrim_pred_postprobe_probDiscrim = lme_function(model_str, postprobe_problemDiscrim, last_day_discrimination, subj_num_mtx, stage_num_mtx);
    last_day_discrim_pred_postprobe_probDiscrim_against = lme_function(model_str_against, postprobe_problemDiscrim, last_day_discrimination, subj_num_mtx, stage_num_mtx)

% does prior probe discrimination predict first-day problem discrimination?
%

    % logistic
    %{
    prior_probe_outside = probe_discr_next_prob_outside; % pre probe only
    first_day_discrimination = first_day_d;
    title_string = [remove_underscore(training_group) ': Prior probe outside predicts first day performance'];
    %plot2var_corr(prior_probe_outside, first_day_discrimination, title_string)
    plot2var(prior_probe_outside, first_day_discrimination, title_string)
    %}

    % normal
    %
    preprobe_problemDiscrim = diff_waits_mean_problem_preprobe; % pre probe only
    first_day_discrimination = first_day_d;
    title_string = [remove_underscore(training_group) ': Prior probe problem discrim predicts first day performance'];
    %plot2var_corr(preprobe_problemDiscrim, first_day_discrimination, title_string)
    plot2var(preprobe_problemDiscrim, first_day_discrimination, title_string)
    tone_discrim_full_pred_firstDay = lme_function(model_str, first_day_discrimination, preprobe_problemDiscrim, subj_num_mtx, stage_num_mtx);
    tone_discrim_full_pred_firstDay_against = lme_function(model_str_against, first_day_discrimination, preprobe_problemDiscrim, subj_num_mtx, stage_num_mtx)
    
% does prior probe discrimination predict last-day problem discrimination?

    % logistic
    %{
    prior_probe_outside = probe_discr_next_prob_outside; % pre probe only
    last_day_discrimination = last_day_d;
    title_string = [remove_underscore(training_group) ': Prior probe outside discrim predicts last day performance'];
    %plot2var_corr(prior_probe_outside, last_day_discrimination, title_string)
    plot2var(prior_probe_outside, last_day_discrimination, title_string)
    %}
    % normal
    %
    preprobe_problemDiscrim = diff_waits_mean_problem_preprobe; % pre probe only
    last_day_discrimination = last_day_d;
    title_string = [remove_underscore(training_group) ': Prior probe problem discrim predicts last day performance'];
    %plot2var_corr(preprobe_problemDiscrim, last_day_discrimination, title_string)
    plot2var(preprobe_problemDiscrim, last_day_discrimination, title_string)
    tone_discrim_full_pred_lastDay = lme_function(model_str, last_day_discrimination, preprobe_problemDiscrim, subj_num_mtx, stage_num_mtx);
    tone_discrim_full_pred_lastDay_against = lme_function(model_str_against, last_day_discrimination, preprobe_problemDiscrim, subj_num_mtx, stage_num_mtx)


% does prior probe discrimination predict number of trials to criterion?

    % logistic
    %{
    prior_probe_outside = probe_discr_past_prob_outside; % pre probe only
    days_to_criterion = days_to_crit;
    title_string = [remove_underscore(training_group) ': Prior probe outside discrim predicts days to criterion'];
    %plot2var_corr(prior_probe_outside, days_to_criterion, title_string)
    plot2var(prior_probe_outside, days_to_criterion, title_string)
    %}

    % normal
    %
    preprobe_problemDiscrim = diff_waits_mean_problem_preprobe; % pre probe only
    days_to_criterion = days_to_crit;
    title_string = [remove_underscore(training_group) ': Prior probe problem discrim predicts days to criterion'];
    %plot2var_corr(preprobe_problemDiscrim, days_to_criterion, title_string)
    plot2var(preprobe_problemDiscrim, days_to_criterion, title_string)   
    tone_discrim_full_pred_days2crit = lme_function(model_str, days_to_criterion, preprobe_problemDiscrim, subj_num_mtx, stage_num_mtx);
    tone_discrim_full_pred_days2crit_against = lme_function(model_str_against, days_to_criterion, preprobe_problemDiscrim, subj_num_mtx, stage_num_mtx)
         
    
% does prior probe problem discrim predict early learning rate
    figure; hold on
    plot(preprobe_problemDiscrim(:), TwoDay_learn_mtx(:), 'o')
    
    window_size = 2.5;
    xpos_probDiscrim = -9:.025:9;
    ypos_probLearn_mean = nan(1,length(xpos_probDiscrim));
    ypos_probLearn_std = nan(1,length(xpos_probDiscrim));
    for iwdw = 1:length(xpos_probDiscrim)
        window_dp = TwoDay_learn_mtx(preprobe_problemDiscrim>xpos_probDiscrim(iwdw)-window_size/2 & preprobe_problemDiscrim<=xpos_probDiscrim(iwdw)+window_size/2);
        ypos_probLearn_mean(iwdw) = nanmean(window_dp);
        ypos_probLearn_std(iwdw) = nanstd(window_dp)./sqrt(sum(~isnan(window_dp)));
    end
    
    plot(xpos_probDiscrim, ypos_probLearn_mean, 'linewidth', 2)
    plot(xpos_probDiscrim, ypos_probLearn_mean-ypos_probLearn_std, 'linewidth', 1)
    plot(xpos_probDiscrim, ypos_probLearn_mean+ypos_probLearn_std, 'linewidth', 1)
    plot(xlim, [1 1].*0, 'k--')
    plot(0.*[1 1], ylim, 'k--')
    
end         
         
%% Figure functions
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
function plot2var_corr(var1, var2, title_str)
    figure
    for iprob = 1:size(var1,2)
        subplot(1,size(var1,2),iprob) 
        hold on
        [r,p] = fit_line(var1(:,iprob), var2(:,iprob)); 
        axis([-10 10 -10 10])
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


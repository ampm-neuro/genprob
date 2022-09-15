function [diff_waits_problem, ratio_waits_problem] = ALL_probe_vs_behavior_terpolate(training_group)
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


%% all probe wait times

% wait times from every probe for every subject (probe,tone,subj)
[probe_wait_times] = ALL_probe_wait_times(training_group, 1); % model fit, interp & extrap, no smooth

% number of probes
num_probes = size(probe_wait_times,1);
num_tones = size(probe_wait_times,2);
num_subjects = size(probe_wait_times,3);


%% mixed model prep
subj_num_mtx = repmat((1:num_subjects)', 1, num_problems); % subject matrix
stage_num_mtx = repmat(1:6, num_subjects, 1); % stage matrix
model_str = 'var1~var2+(1|subject)+(1|problem)';
model_str_against = 'var1~var2+problem+(1|subject)';
model_str_against_double = 'var1~var2+var3+problem+(1|subject)';



%% Use wait times to estimate discrimination of every problem tone pair during every probe test

%
% preallocate (probe, problem, subj)
rich_waits_problem = nan(num_probes, num_problems, num_subjects);
poor_waits_problem = nan(num_probes, num_problems, num_subjects);
diff_waits_problem = nan(num_probes, num_problems, num_subjects);
ratio_waits_problem = nan(num_probes, num_problems, num_subjects);


% iterate through all probes
for isubj = 1:size(probe_wait_times,3)
    for iproblem = 1:num_problems
        for iprobe = 1:num_probes
            
            % tone numbers
            rich_tone_num = all_prob_tone_nums(iproblem, 1);
            poor_tone_num = all_prob_tone_nums(iproblem, 2);
            
            
            % tone wait times (during probe) FULL MODEL
            rich_waits_problem(iprobe, iproblem, isubj) = probe_wait_times(iprobe, rich_tone_num, isubj);
            poor_waits_problem(iprobe, iproblem, isubj) = probe_wait_times(iprobe, poor_tone_num, isubj);
            
            % difference between rich and poor tone wait times (during probe)
            diff_waits_problem(iprobe, iproblem, isubj)...
                = rich_waits_problem(iprobe, iproblem, isubj) - poor_waits_problem(iprobe, iproblem, isubj);
            % ratio between rich and poor tone wait times (during probe)
            ratio_waits_problem(iprobe, iproblem, isubj)...
                = (rich_waits_problem(iprobe, iproblem, isubj) - poor_waits_problem(iprobe, iproblem, isubj))/(rich_waits_problem(iprobe, iproblem, isubj) + poor_waits_problem(iprobe, iproblem, isubj));

        end
    end
end



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






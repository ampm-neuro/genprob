function [pval] = thirds_test(trl_mtx)
% splits trials into highest third and lowest third in terms of tones 
% heard and then compares the thirds with an independent ttest2

% probe trials only
trl_mtx = trl_mtx(trl_mtx(:,3)==0,:);

% number of trials
num_trials = size(trl_mtx,1);
one_third_num_trials = floor(num_trials/3);
two_third_num_trials = ceil((num_trials/3)*2);

% trial tones
trl_tones = floor(trl_mtx(:,2));

% sort tones by frequency
[~, trial_tones_sort_idx] = sort(trl_tones);

% thirds 
low_third_trials = trial_tones_sort_idx(1:one_third_num_trials);
mid_third_trials = trial_tones_sort_idx(one_third_num_trials+1 : two_third_num_trials-1);
high_third_trials = trial_tones_sort_idx(two_third_num_trials:end);

% wait times
low_third_waits = trl_mtx(low_third_trials, 12);
mid_third_waits = trl_mtx(mid_third_trials, 12);
high_third_waits = trl_mtx(high_third_trials, 12);
wait_times = [low_third_waits; mid_third_waits; high_third_waits];

% group labels
group_labels = [ones(size(low_third_waits)); 2.*ones(size(mid_third_waits)); 3.*ones(size(high_third_waits))];

% ttest
[~, pval, ~, stats] = ttest2(low_third_waits, high_third_waits);
tstat = stats.tstat;

% anova
%[pval] = anova1(wait_times, group_labels, 'off');
%fstat = stats.tstat;

%figure; wait_times_plot(trl_mtx);
%figure; errorbar_plot([{low_third_waits}, {mid_third_waits}, {high_third_waits}]);




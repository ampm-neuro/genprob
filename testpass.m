function pass_fail = testpass(trl_mtx, medass_cell)
% 1 for pass 0 for fail


% compute wait durations
[wait_durations_all, wda_freq] = wait_times_prep(trl_mtx,2);

% compute reward probabilities
[prob_dist, pd_freq] = rwd_prob_by_freq(medass_cell);
rich_freq = pd_freq(prob_dist>0.5);
poor_freqs = pd_freq(prob_dist<0.5);

% test differences between rich and poor tone wait times
%

[~, ~, anova_grp] = unique(wda_freq);
anova_p = anovan(wait_durations_all, anova_grp, 'display', 'off');

% posthoc ttests
[~, pval1, ~, stats1] = ttest2(wait_durations_all(wda_freq==poor_freqs(1)), wait_durations_all(wda_freq==rich_freq));
[~, pval2, ~, stats2] = ttest2(wait_durations_all(wda_freq==poor_freqs(2)), wait_durations_all(wda_freq==rich_freq));
tstat = [stats1.tstat stats2.tstat];

% test
if tstat(1)<0 && tstat(2)<0 && pval1<0.05 && pval2<0.05 && anova_p<0.05 && length(wait_durations_all(wda_freq==rich_freq))>=2
    pass_fail = 1;
else
    pass_fail = 0;
end
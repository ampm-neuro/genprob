function [zdiff] = rwd_recency_dep(trl_mtx, nback)
% correlates the proportion of recent trials (in the window spanning nback 
% trials) that were rewarded with the wait time on the current probe

%preallocate
nback_rwd_props = [];
probe_trial_waits = [];
count = 0;

%iterate through trials
for itrl = max(nback)+1:size(trl_mtx)
    
    if trl_mtx(itrl,3)==0
        count = count+1;
        
        % nback trial
        nback_trl = itrl-nback;

        % nback reward history (1 rwd, 0 no rwd)
        nback_rwd_props = [nback_rwd_props; double(any(~isnan(trl_mtx(nback_trl,11))))];
        
        % current wait time
        probe_trial_waits = [probe_trial_waits; trl_mtx(itrl, 12)];

    end
end
prev_rwd_trials = probe_trial_waits(nback_rwd_props==1);
prev_nonrwd_trials = probe_trial_waits(nback_rwd_props==0);

zdiff = (mean(prev_rwd_trials)-mean(prev_nonrwd_trials))/mean([std(prev_rwd_trials) std(prev_nonrwd_trials)]);

%[~,p,~, tstats] = ttest2(probe_trial_waits(nback_rwd_props==1),probe_trial_waits(nback_rwd_props==0));
%t = tstats.tstat;

function [wait_times, frequencies] = wait_times_prep_z(trl_mtx, mean_or_all)
% outputs wait times and corresponding tone frequencies
% for means mean_or_all == 1, all == 2

probe_trials_idx = trl_mtx(:,3)==0;
all_wait_times = trl_mtx(probe_trials_idx,12)+2; %2s for fixed delay
all_tones = trl_mtx(probe_trials_idx,2);

% zscore waits
all_wait_times = zscore_mtx(all_wait_times);

if mean_or_all == 1 % means onlynbnb
    frequencies = unique(all_tones);
    wait_times = nan(size(frequencies));
    for itone = 1:length(frequencies)
        tone = frequencies(itone);
        wait_times(itone) = mean(all_wait_times(all_tones==tone));
    end
else % all waits
    wait_times = all_wait_times;
    frequencies = all_tones;
end
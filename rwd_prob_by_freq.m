function [rwd_prob, unq_frq] = rwd_prob_by_freq(medass_cell)
% outputs the probability of reward for each frequency

% raw frequencies
freq_raw = floor(medass_cell{21});

% raw reward probablities
rwd_probs_raw = medass_cell{22};

% unique freq
unq_frq = unique(freq_raw);

% compute reward probabilities
rwd_prob = nan(size(unq_frq));
for ifrq = 1:length(unq_frq)
    rwd_prob(ifrq) = mean(rwd_probs_raw(freq_raw==unq_frq(ifrq)));
end

% ensure proportion (not probability)
if rwd_prob(1)>1
    rwd_prob = rwd_prob./100;
end

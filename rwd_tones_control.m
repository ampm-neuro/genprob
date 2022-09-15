% script to generate the tones used for training

%load('10trial_freq_idx_3tones_ctrl')
load('rwd_tones_7_test.mat')

%compute number of tones indicated by tone_stimuli_heatmap
num_tones = length([unqfrq41(tsh_7_ctl(1,:)==2)...
            unqfrq41(tsh_7_ctl(1,:)==1)]);
num_trials = size(tsh_7_ctl,1);

%preallocate training frequency matrix
%rwd1_nrwd_seqence = nan(num_trials,num_tones);
rwd1_nrwd_seqence = nan(num_trials,7);

%load training frequencies
for i = 1:num_trials
    
  tones = [unqfrq41(tsh_7_ctl(i,:)==2)...
            unqfrq41(tsh_7_ctl(i,:)==1)];
        
  rwd1_nrwd_seqence(i,1:length(tones)) = [unqfrq41(tsh_7_ctl(i,:)==2)...
            unqfrq41(tsh_7_ctl(i,:)==1)];
        
end

[bar_vect_non_rwd_trls, all_bar_vect] = rwd_tone_seq_sim(rwd1_nrwd_seqence, unqfrq41, probe_hm2);

clearvars i ans
all_bar_vect = [all_bar_vect(1:length(all_bar_vect)/2); all_bar_vect(length(all_bar_vect)/2 + 1 : end)];
all_bar_vect = sum(all_bar_vect);
tone_idx = []; 
for itone = 1:length(all_bar_vect) 
    if all_bar_vect(itone)>0 
        tone_idx = [tone_idx; repmat(itone,all_bar_vect(itone),1)]; 
    end
end

new_rwd_tones = tone_idx(2:3:end);

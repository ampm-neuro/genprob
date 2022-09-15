% script to generate the tones used for training

%load('10trial_freq_idx_3tones')
load('rwd_tones_7_test.mat')


%compute number of tones indicated by tone_stimuli_heatmap
num_tones = length([unqfrq41(tsh_7(1,:)==2)...
            unqfrq41(tsh_7(1,:)==1)]);
num_trials = size(tsh_7,1);

%preallocate training frequency matrix
%rwd1_nrwd_seqence = nan(num_trials,num_tones);
rwd1_nrwd_seqence = nan(num_trials,10);

%load training frequencies
for i = 1:num_trials
    
  tones = [unqfrq41(tsh_7(i,:)==2)...
            unqfrq41(tsh_7(i,:)==1)];
        
  rwd1_nrwd_seqence(i,1:length(tones)) = [unqfrq41(tsh_7(i,:)==2)...
            unqfrq41(tsh_7(i,:)==1)];
        
end

rwd_tone_seq_sim(rwd1_nrwd_seqence, unqfrq41);

clearvars i ans

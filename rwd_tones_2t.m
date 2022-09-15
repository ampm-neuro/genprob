% script to generate the tones used for training

load('rich_tones_mtx_2.mat')
load('unqfrq41.mat')

%compute number of tones indicated by tone_stimuli_heatmap
num_tones = length([unqfrq41(rich_tones_mtx(1,:)==2)...
            unqfrq41(rich_tones_mtx(1,:)==1)]);
num_trials = size(rich_tones_mtx,1);

%preallocate training frequency matrix
rwd1_nrwd_seqence = nan(num_trials,2);
rwd1_nrwd_seqence_ctl = nan(num_trials,2);

%load training frequencies
for i = 1:num_trials
    
  tones = [unqfrq41(rich_tones_mtx(i,:)==2)...
            unqfrq41(rich_tones_mtx(i,:)==1)];
        
  rwd1_nrwd_seqence(i,1:length(tones)) = [unqfrq41(rich_tones_mtx(i,:)==2)...
            unqfrq41(rich_tones_mtx(i,:)==1)];
        
end
rwd_tone_seq_sim(rwd1_nrwd_seqence, unqfrq41);


%load control (scrambled) frequencies
for i = 1:num_trials
    
  tones = [unqfrq41(rich_tones_mtx_ctl(i,:)==2)...
            unqfrq41(rich_tones_mtx_ctl(i,:)==1)];
        
  rwd1_nrwd_seqence_ctl(i,1:length(tones)) = [unqfrq41(rich_tones_mtx_ctl(i,:)==2)...
            unqfrq41(rich_tones_mtx_ctl(i,:)==1)];
        
end
rwd_tone_seq_sim(rwd1_nrwd_seqence_ctl, unqfrq41);

clearvars i ans

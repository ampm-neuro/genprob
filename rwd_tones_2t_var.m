% script to generate the tones used for training

load('rich_tones_mtx_2_var_5.mat')
load('unqfrq41.mat')

%compute number of tones indicated by tone_stimuli_heatmap
num_tones = length([unqfrq41(rich_tones_mtx_novar(1,:)==2)...
            unqfrq41(rich_tones_mtx_novar(1,:)==1)]);
num_trials = size(rich_tones_mtx_novar,1);

%preallocate training frequency matrix
rwd1_nrwd_seqence_novar = nan(num_trials,2);
rwd1_nrwd_seqence_lovar = nan(num_trials,2);
rwd1_nrwd_seqence_mevar = nan(num_trials,2);
rwd1_nrwd_seqence_hivar = nan(num_trials,2);
rwd1_nrwd_seqence_exvar = nan(num_trials,2);

%load training frequencies NO VARIANCE
for i = 1:num_trials
    
  tones = [unqfrq41(rich_tones_mtx_novar(i,:)==2)...
            unqfrq41(rich_tones_mtx_novar(i,:)==1)];
        
  rwd1_nrwd_seqence_novar(i,1:length(tones)) = [unqfrq41(rich_tones_mtx_novar(i,:)==2)...
            unqfrq41(rich_tones_mtx_novar(i,:)==1)];
        
end
rwd_tone_seq_sim(rwd1_nrwd_seqence_novar, unqfrq41);


%load control frequencies LOW VARIANCE
for i = 1:num_trials
    
  tones = [unqfrq41(rich_tones_mtx_lovar(i,:)==2)...
            unqfrq41(rich_tones_mtx_lovar(i,:)==1)];
        
  rwd1_nrwd_seqence_lovar(i,1:length(tones)) = [unqfrq41(rich_tones_mtx_lovar(i,:)==2)...
            unqfrq41(rich_tones_mtx_lovar(i,:)==1)];
        
end
rwd_tone_seq_sim(rwd1_nrwd_seqence_lovar, unqfrq41);


%load control frequencies MEDIUM VARIANCE
for i = 1:num_trials
    
  tones = [unqfrq41(rich_tones_mtx_mevar(i,:)==2)...
            unqfrq41(rich_tones_mtx_mevar(i,:)==1)];
        
  rwd1_nrwd_seqence_mevar(i,1:length(tones)) = [unqfrq41(rich_tones_mtx_mevar(i,:)==2)...
            unqfrq41(rich_tones_mtx_mevar(i,:)==1)];
        
end
rwd_tone_seq_sim(rwd1_nrwd_seqence_mevar, unqfrq41);


%load control frequencies HIGH VARIANCE
for i = 1:num_trials
    
  tones = [unqfrq41(rich_tones_mtx_hivar(i,:)==2)...
            unqfrq41(rich_tones_mtx_hivar(i,:)==1)];
        
  rwd1_nrwd_seqence_hivar(i,1:length(tones)) = [unqfrq41(rich_tones_mtx_hivar(i,:)==2)...
            unqfrq41(rich_tones_mtx_hivar(i,:)==1)];
        
end
rwd_tone_seq_sim(rwd1_nrwd_seqence_hivar, unqfrq41);


%load control frequencies EXTREME VARIANCE
for i = 1:num_trials
    
  tones = [unqfrq41(rich_tones_mtx_exvar(i,:)==2)...
            unqfrq41(rich_tones_mtx_exvar(i,:)==1)];
        
  rwd1_nrwd_seqence_exvar(i,1:length(tones)) = [unqfrq41(rich_tones_mtx_exvar(i,:)==2)...
            unqfrq41(rich_tones_mtx_exvar(i,:)==1)];
        
end
rwd_tone_seq_sim(rwd1_nrwd_seqence_exvar, unqfrq41);

clearvars i ans

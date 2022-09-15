rwd_tones_2t

load('unqfrq41_gain_Y2019_M08_D27.mat')

numtones = [10 10];

% training
for train = 1:size(rwd1_nrwd_seqence,1) 
    tonecsv = [];  
    
    for itone = 1:size(rwd1_nrwd_seqence,2)
        tonecsv = [tonecsv repmat(rwd1_nrwd_seqence(train,itone), [1, numtones(itone)])] ;
    end
    train_str = num2str(train);
    if length(train_str)<2; train_str = ['0' train_str]; end
    
    csvwrite(['tones_train' num2str(train) '.csv'], tonecsv)
    
    for ivol = 1:length(tonecsv)
        tonecsv_vol(ivol) = vols(unqfrq41(1,:)==tonecsv(ivol));
    end
    tonecsv_vol = round(tonecsv_vol);
    
    csvwrite(['vol_train' num2str(train) '.csv'], tonecsv_vol)
    
end


% control
for train = 1:size(rwd1_nrwd_seqence_ctl,1) 
    tonecsv = [];  
    
    for itone = 1:size(rwd1_nrwd_seqence_ctl,2)
        tonecsv = [tonecsv repmat(rwd1_nrwd_seqence_ctl(train,itone), [1, numtones(itone)])]; 
    end
    train_str = num2str(train);
    if length(train_str)<2; train_str = ['0' train_str]; end
    
    csvwrite(['tones_ctrl' num2str(train) '.csv'], tonecsv)
    
    for ivol = 1:length(tonecsv)
        tonecsv_vol(ivol) = vols(unqfrq41(1,:)==tonecsv(ivol));
    end
    tonecsv_vol = round(tonecsv_vol);
    
    csvwrite(['vol_ctrl' num2str(train) '.csv'], tonecsv_vol)
    
end
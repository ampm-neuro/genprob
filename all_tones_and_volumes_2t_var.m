rwd_tones_2t_var2

load('unqfrq41_gain_Y2019_M11_D16.mat')

numtones = [10 10];

% no var
%{
for train = 1:size(rwd1_nrwd_seqence_novar,1) 
    tonecsv = [];  
    
    for itone = 1:size(rwd1_nrwd_seqence_novar,2)
        tonecsv = [tonecsv repmat(rwd1_nrwd_seqence_novar(train,itone), [1, numtones(itone)])] ;
    end
    train_str = num2str(train);
    if length(train_str)<2; train_str = ['0' train_str]; end
    
    csvwrite(['tones_novar' num2str(train) '.csv'], tonecsv)
    
    for ivol = 1:length(tonecsv)
        tonecsv_vol(ivol) = vols(unqfrq41(1,:)==tonecsv(ivol));
    end
    tonecsv_vol = round(tonecsv_vol);
    
    csvwrite(['vol_novar' num2str(train) '.csv'], tonecsv_vol)
    
end
%}


% low var
for train = 1:size(rwd1_nrwd_seqence_lovar,1) 
    tonecsv = [];  
    
    for itone = 1:size(rwd1_nrwd_seqence_lovar,2)
        tonecsv = [tonecsv repmat(rwd1_nrwd_seqence_lovar(train,itone), [1, numtones(itone)])]; 
    end
    train_str = num2str(train);
    if length(train_str)<2; train_str = ['0' train_str]; end
    
    csvwrite(['tones_lovar' num2str(train) '.csv'], tonecsv)
    
    for ivol = 1:length(tonecsv)
        tonecsv_vol(ivol) = vols(unqfrq41(1,:)==tonecsv(ivol));
    end
    tonecsv_vol = round(tonecsv_vol);
    
    csvwrite(['vol_lovar' num2str(train) '.csv'], tonecsv_vol)
    
end


% medium var
for train = 1:size(rwd1_nrwd_seqence_mevar,1) 
    tonecsv = [];  
    
    for itone = 1:size(rwd1_nrwd_seqence_mevar,2)
        tonecsv = [tonecsv repmat(rwd1_nrwd_seqence_mevar(train,itone), [1, numtones(itone)])]; 
    end
    train_str = num2str(train);
    if length(train_str)<2; train_str = ['0' train_str]; end
    
    csvwrite(['tones_mevar' num2str(train) '.csv'], tonecsv)
    
    for ivol = 1:length(tonecsv)
        tonecsv_vol(ivol) = vols(unqfrq41(1,:)==tonecsv(ivol));
    end
    tonecsv_vol = round(tonecsv_vol);
    
    csvwrite(['vol_mevar' num2str(train) '.csv'], tonecsv_vol)
    
end


% high var
for train = 1:size(rwd1_nrwd_seqence_hivar,1) 
    tonecsv = [];  
    
    for itone = 1:size(rwd1_nrwd_seqence_hivar,2)
        tonecsv = [tonecsv repmat(rwd1_nrwd_seqence_hivar(train,itone), [1, numtones(itone)])]; 
    end
    train_str = num2str(train);
    if length(train_str)<2; train_str = ['0' train_str]; end
    
    csvwrite(['tones_hivar' num2str(train) '.csv'], tonecsv)
    
    for ivol = 1:length(tonecsv)
        tonecsv_vol(ivol) = vols(unqfrq41(1,:)==tonecsv(ivol));
    end
    tonecsv_vol = round(tonecsv_vol);
    
    csvwrite(['vol_hivar' num2str(train) '.csv'], tonecsv_vol)
    
end

% ex var
%{
for train = 1:size(rwd1_nrwd_seqence_exvar,1) 
    tonecsv = [];  
    
    for itone = 1:size(rwd1_nrwd_seqence_exvar,2)
        tonecsv = [tonecsv repmat(rwd1_nrwd_seqence_exvar(train,itone), [1, numtones(itone)])]; 
    end
    train_str = num2str(train);
    if length(train_str)<2; train_str = ['0' train_str]; end
    
    csvwrite(['tones_exvar' num2str(train) '.csv'], tonecsv)
    
    for ivol = 1:length(tonecsv)
        tonecsv_vol(ivol) = vols(unqfrq41(1,:)==tonecsv(ivol));
    end
    tonecsv_vol = round(tonecsv_vol);
    
    csvwrite(['vol_exvar' num2str(train) '.csv'], tonecsv_vol)
    
end
%}
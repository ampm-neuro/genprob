
subject = '699438m1';

%
for iprob = 7:9
    for iday = 1:2

        vp = ['F:\ampm\data\' subject '\gen0' num2str(iprob) '\d0' num2str(iday)]
        [all_pkl_frames, trl_idx] = process_splitvids(vp); 
        title(['problem ' num2str(iprob) '; day ' num2str(iday)]); 
        drawnow
    end
end
%}
%
for iprob = 10:12
    for iday = 1:2
        vp = ['F:\ampm\data\' subject '\gen' num2str(iprob) '\d0' num2str(iday)] 
        [all_pkl_frames, trl_idx] = process_splitvids(vp); 
        title(['problem ' num2str(iprob) '; day ' num2str(iday)]); 
        drawnow
    end
end
%}

%
for iprobe = 1:8
    iprobe
    vp = ['F:\ampm\data\' subject '\probe0' num2str(iprobe) '\d01'] 
    [all_pkl_frames, trl_idx] = process_splitvids(vp); 
    title(['probe ' num2str(iprobe)]); 
    drawnow
end
%}
%
for NT = 1:2
    try
    vp = ['F:\ampm\data\' subject '\NT0' num2str(NT) '\d01']; 
    [all_pkl_frames, trl_idx] = process_splitvids(vp); 
    title(['NT 0' num2str(NT)]); 
    drawnow
    catch
    end
end
%}


for iprob = 12
    for iday = 1
        vp = ['F:\ampm\data\' subject '\gen' num2str(iprob) '\d0' num2str(iday)] 
        [all_pkl_frames, trl_idx] = process_splitvids(vp); 
        title(['problem ' num2str(iprob) '; day ' num2str(iday)]); 
        drawnow
    end
end


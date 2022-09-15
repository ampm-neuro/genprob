function [reliability, trial_mtx]  = one_cell_reliability(neurons, session_mtx_cell)
% load session_mtx_cell computed inside image_pop_mds 

% preallocate
reliability = nan(length(session_mtx_cell), 1);
trial_mtx = nan(41, 41, length(session_mtx_cell));

% iterate through neurons
for cell_num = neurons
    
    for isesh = 1:length(session_mtx_cell)
        activity = session_mtx_cell{isesh}{cell_num};
        
        % limit trials
        activity = activity(1:41,:);
        
        act_rel = nan(size(activity,1));
        for itrl1 = 1:size(activity,1)
            for itrl2 = 1:size(activity,1)
                if itrl1>=itrl2
                    continue
                end
                act_rel(itrl1,itrl2) = corr(activity(itrl1,:)', activity(itrl2,:)');
            end
        end
        
        reliability(isesh) = nanmean( act_rel(:) );
        
        % load if trials are limited to 41
        trial_mtx(:,:,isesh) = act_rel;
        
    end
    
    
end
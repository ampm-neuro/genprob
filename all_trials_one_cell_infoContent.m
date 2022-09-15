function info_content = all_trials_one_cell_infoContent(neurons, session_mtx_cell)
% load session_mtx_cell computed inside image_pop_mds 

% number of trials in each session
session_trials = nan(length(session_mtx_cell),1);
for isesh = 1:length(session_mtx_cell)
    session_trials(isesh) = size(session_mtx_cell{isesh}{1},1);
end
cum_session_trials = cumsum(session_trials);

% iterate through neurons
for cell_num = neurons
    fullcell = []; 
    
    for isesh = 1:length(session_mtx_cell)
        %fullcell = [fullcell; session_mtx_cell{isesh}{cell_num}(1:41,:)];
        fullcell = [fullcell; session_mtx_cell{isesh}{cell_num}];
    end
    

    % norm rates
    fullcell = norm_mtx(fullcell')';
    
    %figure; hold on
    
    
    % compute info content
    info_content = nan(size(fullcell,1), 1);
    for itrl = 1:size(fullcell,1)
        if sum(fullcell(itrl),2) ~= 0 && sum(~isnan(fullcell(itrl)),2) ~= 0

            info_content(itrl) = contrast_score(fullcell(itrl,:));

        end
    end
    
    %info_content
    %figure; plot(info_content)
    %hold on; plot(norm_mtx(info_content))
    
end
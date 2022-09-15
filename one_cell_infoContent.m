function info_content = one_cell_infoContent(neuron, session_mtx_cell)
% load session_mtx_cell computed inside image_pop_mds 

% preallocate
info_content = nan(1, length(session_mtx_cell));

for isesh = 1:length(session_mtx_cell)

    % get activity for this cell on this trial
    trials_activity = session_mtx_cell{isesh}{neuron}(1:41,:);


    % normalize activity
    trials_activity = norm_mtx(trials_activity')';

    % compute info content
    info_content_trials = nan(size(trials_activity,1),1);
    for itrl = 1:size(trials_activity,1)

        % if there is activity and it's not all nans
        if sum(trials_activity(itrl,:),2) ~= 0 && sum(~isnan(trials_activity(itrl,:)),2) ~= 0
            info_content_trials(itrl) = contrast_score(trials_activity(itrl,:));
        end
    end

    % load average across all trials
    info_content(isesh) = mean(info_content_trials);
    
end
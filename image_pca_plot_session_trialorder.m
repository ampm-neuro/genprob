function image_pca_plot_session_trialorder(session_paths)
% plot all time bins colored by trial


if ~iscell(session_paths)
    session_paths_hold{1} = session_paths;
    session_paths = session_paths_hold;
end


tses = [0.2 1.1 1.0 2.0 8.5 2.0];

% all_matrices is cell X time X trial



%% compute activity matrix
all_matrices = [];
for isesh = 1:length(session_paths)
    trial_activity_mtx = image_trial_activity_mtx(session_paths{isesh});
    if size(trial_activity_mtx,3)>size(all_matrices,3) && ~isempty(all_matrices)
        all_matrices = cat(1, all_matrices, trial_activity_mtx(:,:,1:size(all_matrices,3)));
    elseif size(trial_activity_mtx,3)>size(all_matrices,3) && ~isempty(all_matrices)
        all_matrices = cat(1, all_matrices(:,:,size(trial_activity_mtx,3)), trial_activity_mtx);
    else
        all_matrices = cat(1, all_matrices, trial_activity_mtx);
    end
end
all_matrices = all_matrices(:,~isnan(sum(all_matrices(:,:,1))),:);



%% prepare input for dim red

% remove empty time bins and empty trials
%tones = tones(sum(sum(~isnan(all_matrices)),2)>0);
all_matrices = all_matrices(:,sum(sum(~isnan(all_matrices)),3)>0, sum(sum(~isnan(all_matrices)),2)>0);

% remove time windows after end of random delay
event_frame = cumsum(tses).*100;
all_matrices = all_matrices(:, 1:event_frame(4) ,:);
size(all_matrices)

% zscore each cell
%
for ic = 1:size(all_matrices,1)
    hold_ic = all_matrices(ic,:,:);
    all_matrices(ic,:,:) = reshape(zscore_mtx(hold_ic(:)), [1, size(all_matrices,2), size(all_matrices,3)]);
end
%}



%% dim reduction
% be sure to add tensor toolbox to path
am_model = cp_als(tensor(all_matrices), 10); % fit CP model with 10 components
am_model = double(am_model);



%% stage index to 1:4 (Nose Poke, Movement, Head Entry, Tone On)
stage_idx = 1:size(am_model,2);
event_frame = cumsum(tses).*100;
for istage = 1:4
    stage_idx(stage_idx>=istage & stage_idx<=event_frame(istage)) = istage;
end




%% plot

dims = [1 2 3];
colors = parula(size(am_model,3));

%
for istage = 1:4
    figure; hold on; 
    
    for itrl=1:size(am_model,3)

        % indices
        istage_idx = stage_idx==istage;
        
        % plot transparent
        scatter1 = scatter3(am_model(dims(1), ~istage_idx, itrl), am_model(dims(2), ~istage_idx, itrl), am_model(dims(3), ~istage_idx, itrl), 8, colors(itrl,:), 'filled');  
        alpha(scatter1, 0.05)
        
        % plot stage
        scatter3(am_model(dims(1), istage_idx, itrl), am_model(dims(2), istage_idx, itrl), am_model(dims(3), istage_idx, itrl), 8, colors(itrl,:), 'filled'); 

    end
    
    grid on
    axis square
    title(['Stage ' num2str(istage)])
    set(gca,'TickLength',[0, 0]); box off;
end
%}























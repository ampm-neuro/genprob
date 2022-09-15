function [all_matrices] = image_pca_plot_session_trialorder_unqcells(subject_ids, mevar_or_hivar)
% plot all time bins colored by trial



tses = [0.2 1.1 1.0 2.0 8.5 2.0];



%% compute activity matrix

all_matrices = [];
for isubj = 1:length(subject_ids)
    
    % cell_regist_mtx
    clearvars cell_regist_mtx
    load(['cell_reg_' subject_ids{isubj} '.mat']); cell_regist_mtx = cell_registered_struct.cell_to_index_map; 
    
    % chron
    scri = session_chron_reorder_idx;
    cell_regist_mtx = cell_regist_mtx(:,scri);
    
    % get session paths
    last_prob_sessions = [];
    for iprob = 1:6
        prob_sessions = get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_' mevar_or_hivar '_imaging_hpc\' subject_ids{isubj} ], ['var0' num2str(iprob)], 'LED');
        last_prob_sessions = [last_prob_sessions; prob_sessions(end)];
    end
    session_cell = [...
        get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_' mevar_or_hivar '_imaging_hpc\' subject_ids{isubj} ], 'preprobe', 'LED');...
        get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_' mevar_or_hivar '_imaging_hpc\' subject_ids{isubj} ] , 'postprobe', 'LED');...
        get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_' mevar_or_hivar '_imaging_hpc\' subject_ids{isubj} ], 'var0', '-01', 'LED');...
        last_prob_sessions];
    session_cell = session_cell(scri)
    
    % build super matrix of all activity
    all_matrices_subj = [];
    for isesh = 1:4%length(session_cell)-3 : length(session_cell)%length(session_cell)
        session_cell{isesh}
        
        
        % identify cells whos last appearance is in this session
        %last_appearance_cellnums = cell_regist_mtx(cell_regist_mtx(:,isesh)>0 & sum(cell_regist_mtx(:,isesh:end)>0,2)==1, isesh);
        %if isempty(last_appearance_cellnums); continue; end
        
        % skip nonprobe
        %if ~contains(session_cell{isesh}, 'preprobe_02')% || ~contains(session_cell{isesh}, 'postprobe_02')
        %    continue
        %end
        
        % compute session mtx
        trial_activity_mtx = image_trial_activity_mtx(session_cell{isesh});
                
        % only keep last-appearance cells
        %trial_activity_mtx = trial_activity_mtx(last_appearance_cellnums,:,:);
        
        % combine across sessions
        if size(trial_activity_mtx,3)>size(all_matrices_subj,3) && ~isempty(all_matrices_subj)
            all_matrices_subj = cat(1, all_matrices_subj, trial_activity_mtx(:,:,1:size(all_matrices_subj,3)));
        elseif size(trial_activity_mtx,3)<size(all_matrices_subj,3) && ~isempty(all_matrices_subj)
            all_matrices_subj = cat(1, all_matrices_subj(:,:,1:size(trial_activity_mtx,3)), trial_activity_mtx);
        else
            all_matrices_subj = cat(1, all_matrices_subj, trial_activity_mtx);
        end
    end
    
    
    
    % combine across subjects
    if size(all_matrices_subj,3)>size(all_matrices,3) && ~isempty(all_matrices)
        all_matrices = cat(1, all_matrices, all_matrices_subj(:,:,1:size(all_matrices,3)));
    elseif size(all_matrices_subj,3)<size(all_matrices,3) && ~isempty(all_matrices)
        all_matrices = cat(1, all_matrices(:,:,1:size(all_matrices_subj,3)), all_matrices_subj);
    else
        all_matrices = cat(1, all_matrices, all_matrices_subj);
    end
    
end



%% prepare input for dim red

% remove empty time bins and empty trials
%tones = tones(sum(sum(~isnan(all_matrices)),2)>0);
all_matrices = all_matrices(:, sum(sum(~isnan(all_matrices)),3)>0, sum(sum(~isnan(all_matrices)),2)>0);

% remove time windows after end of random delay
event_frame = cumsum(tses).*100;
all_matrices = all_matrices(:, 1:event_frame(4) ,:);

% zscore each cell
%{
for ic = 1:size(all_matrices,1)
    hold_ic = all_matrices(ic,:,:);
    all_matrices(ic,:,:) = reshape(zscore_mtx(hold_ic(:)), [1, size(all_matrices,2), size(all_matrices,3)]);
end
%}

% zscore each cell EACH TRIAL
%
for ic = 1:size(all_matrices,1)
    for itrl = 1:size(all_matrices,3)
        hold_ic = all_matrices(ic,:,itrl);
        all_matrices(ic,:,itrl) = zscore_mtx(hold_ic(:));
    end
end
%}

% remove cells that don't handle the zscoring procedure
all_matrices = all_matrices(~isnan(sum(sum(all_matrices,3),2)), :, :);

size(all_matrices)


%% dim reduction
% be sure to add tensor toolbox to path
addpath(genpath('C:\Users\ampm1\Documents\MATLAB\generalization\tensor_toolbox-v3.2.1'))
am_model = cp_als(tensor(all_matrices), 10); % fit CP model with 10 components
am_model = double(am_model);



%% stage index to 1:5 (Nose Poke, Movement, Head Entry, Tone On, wait)
stage_idx = 1:size(am_model,2);
event_frame = cumsum(tses).*100;
for istage = 1:4
    stage_idx(stage_idx>=istage & stage_idx<=event_frame(istage)) = istage;
end




%% plot

dims = [1 2 3];
colors = parula(size(am_model,3));

%{
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

% plot overall (no opacity)
figure; hold on; 
for itrl=1:size(am_model,3)
    %scatter3(am_model(dims(1), :, itrl), am_model(dims(2), :, itrl), am_model(dims(3), :, itrl), 8, colors(itrl,:), 'filled');
    plot3(am_model(dims(1), :, itrl), am_model(dims(2), :, itrl), am_model(dims(3), :, itrl), '-', 'color', colors(itrl,:));
end

grid on
axis square
title(['All stages'])
set(gca,'TickLength',[0, 0]); box off;
%}























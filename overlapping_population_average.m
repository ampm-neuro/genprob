function [overlap_mtx, probe_mtx, probe_lines] = overlapping_population_average(subject_ids)

% preallocate
overlap_mtx = [];


% subject folders
cell_reg_data_folder_path = 'C:\Users\ampm1\Documents\MATLAB\generalization\cellreg_data';
session_folders = get_folder_paths_all(cell_reg_data_folder_path, 0);

% preallocate cell counts
cell_cts = nan(length(subject_ids),1);

% iterate through subjects
for isubj = 1:length(subject_ids)
    subj_local = subject_ids{isubj};
    
    % subject folder
    subj_folder = session_folders(contains(session_folders, subj_local));
    subj_folder = subj_folder{1};
    
    % get cell_reg file
    cell_reg_files = get_file_paths_targeted(subj_folder , 'cell_reg_', '.mat');
    cell_reg_file = cell_reg_files{1};
    
    % load
    load(cell_reg_file, 'cell_registered_struct')
    
    % compute matrix
    cell_regist_mtx = cell_registered_struct.cell_to_index_map;
    cell_cts(isubj) = size(cell_regist_mtx,1);
    
    % correct for missing sessions
    if strcmp(subject_ids{isubj}, '690330m1')
        cell_regist_mtx = [cell_regist_mtx(:,1:9) nan(size(cell_regist_mtx(:,1))) cell_regist_mtx(:,10:end)];
    end
    
    session_chron_reorder_idx = [1 9 15 2 10 16 3 11 17 4 12 18 5 13 19 6 14 20 7 8];
    com = cell_overlap_mtx({cell_regist_mtx}, session_chron_reorder_idx, subj_local);
    
    % set diaganal to nan
    com(logical(eye(length(com)))) = nan;
    
    % load
    overlap_mtx = cat(3, overlap_mtx, com);
end

% animal-level mean
figure; imagesc(nanmean(overlap_mtx,3));
axis square;
title mean;
caxis([0 .5])
colorbar

% weighted mean (by number of unique cells)
%{
overlap_mtx_weighted = nan(size(overlap_mtx,1), size(overlap_mtx,2), sum(cell_cts));
loval = 0;
for isubj=1:length(cell_cts)
    overlap_mtx_weighted(:,:,loval+1:loval+cell_cts(isubj)) = repmat(overlap_mtx(:,:,isubj), [1,1,cell_cts(isubj)]);
    loval = cell_cts(isubj);
end
figure; imagesc(nanmean(overlap_mtx_weighted,3));
axis square;
title weightedmean;
caxis([0 .5])
colorbar
%}

% probes only
probe_mtx = nan(8,8, size(overlap_mtx,3));
probenums = [1:3:19 20];
for iprobe1 = 1:8
    for iprobe2 = 1:8
        probe_mtx(iprobe1, iprobe2, :) = overlap_mtx(probenums(iprobe1), probenums(iprobe2), :);
    end
end
figure; imagesc(nanmean(probe_mtx,3));
axis square;
title mean;
caxis([0 .5])
colorbar

% probes errorbarplot
probe_lines = nan(size(probe_mtx,3), 7);
for isubj = 1:size(probe_mtx,3)
    for iprobe1 = 1:7
            probe_lines(isubj, iprobe1) = probe_mtx(iprobe1, iprobe1+1, isubj);
    end
end

figure; hold on; 
errorbar_mtx(probe_lines)






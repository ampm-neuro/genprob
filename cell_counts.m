function [session_ct_mtx, unique_ct_vect] = cell_counts(subject_ids)
% output session cell counts and unique cell counts

% preallocate
session_ct_mtx = nan(length(subject_ids), 20);
unique_ct_vect = nan(length(subject_ids),1);

% organize sessions into chronological order
session_chron_reorder_idx = [1 9 15 2 10 16 3 11 17 4 12 18 5 13 19 6 14 20 7 8];

% subject folders
cell_reg_data_folder_path = 'C:\Users\ampm1\Documents\MATLAB\generalization\cellreg_data';
session_folders = get_folder_paths_all(cell_reg_data_folder_path, 0);

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
    
    % cell map index
    cell_regist_mtx = cell_registered_struct.cell_to_index_map;
    
    % correct for missing sessions
    if strcmp(subject_ids{isubj}, '690330m1')
        cell_regist_mtx = [cell_regist_mtx(:,1:9) nan(size(cell_regist_mtx(:,1))) cell_regist_mtx(:,10:end)];
    end
    
    % sort chron
    cell_regist_mtx = cell_regist_mtx(:, session_chron_reorder_idx);
    
    % session counts
    session_ct_mtx(isubj, :) = sum(cell_regist_mtx>0);
    
    % unique counts
    unique_ct_vect(isubj) = size(cell_regist_mtx,1);
end




end
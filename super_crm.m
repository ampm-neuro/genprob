function [super_crm_out, subject_crm_cell] = super_crm(subj_ids)
% makes a single cell reg matrix used in some plotting functions by opening
% each subjects saved cell_reg file and stacking the crms on top of each
% other (after renumbering to keep 1 number per unique cell)

% preallocate
super_crm_out = [];
subject_crm_cell = cell(1, length(subj_ids));

% iterate through subjects
for isubject = 1:length(subj_ids)
    
    % load cell reg file
    load(['cellreg_data\' subj_ids{isubject} '\cell_reg_' subj_ids{isubject} '.mat'], 'cell_registered_struct');
    
    % cell registration matrix
    if strcmp(subj_ids{isubject}, '690330m1')
        cell_regist_mtx = cell_registered_struct.cell_to_index_map;
        cell_regist_mtx = [cell_regist_mtx(:,1:9) nan(size(cell_regist_mtx(:,1))) cell_regist_mtx(:,10:end)];
    else
        cell_regist_mtx = cell_registered_struct.cell_to_index_map;
    end

    % load    
    super_crm_out = [super_crm_out; cell_regist_mtx];
    subject_crm_cell{isubject} = cell_regist_mtx;

end



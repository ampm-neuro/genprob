function subj_ids = unique_subjs(varargin)
% find unique subjects in behavioral folder that also have files with the
% included strings

   
% behavioral folder
beh_fld = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\';

% all folders
if ~isempty(varargin)
    all_files = get_file_paths_targeted(beh_fld, varargin);
else
    all_files = get_file_paths_targeted(beh_fld);
end

% all subjects
all_subj_ids = [];
for ifile = 1:size(all_files,1)
    all_subj_ids = [all_subj_ids; find_subj_id(all_files{ifile})];
end
    
% unique
subj_ids = unique(all_subj_ids, 'rows');
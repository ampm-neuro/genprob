function [cell_registered_struct, cregstruct_meta, cregstruct_intermediates] = footprints_meta_reverse
%footprints_meta_reverse recreates original 20 session footprint files
% using the corrected files computed with the meta footprints_meta 
% procedure and subsequent cell reg
%
%   the goal of the meta procedure is the use cellreg independently 
%   on groups of sessions, and then use cell reg again on corrected  
%   footprints output by those cell reg procedures 
% 
%   this function takes the final footprints (corrected) and reforms the 
%   original 20 footprint files
%
%   cregstruct_meta is the cell_registered_struct from the final meta cell 
%   reg cregstruct_intermediates is a cell containing the
%   cell_registered_struct from each of the intermediate cell reg
%   procedures. Each cell corresponds to a column in 
%   cregstruct_meta.cell_to_index_map
%
%   if an intermediate stage grouping only had 1 session, leave the
%   corresponding cregstruct_intermediates cell empty
%   
%   
% 

% specific for subject 683472m3
% load
load('C:\Users\ampm1\Desktop\cellreg_temp\683472m3\nested\take1\cellRegistered_meta.mat');
cregstruct_meta = cell_registered_struct;


cregstruct_intermediates = cell(1,4); % 4th cell is empty
load('C:\Users\ampm1\Desktop\cellreg_temp\683472m3\nested\g07d02-p4\cellRegistered_fin.mat');
cregstruct_intermediates{1} = cell_registered_struct;
load('C:\Users\ampm1\Desktop\cellreg_temp\683472m3\nested\g10d02-p8\cellRegistered_fin.mat');
cregstruct_intermediates{2} = cell_registered_struct;
load('C:\Users\ampm1\Desktop\cellreg_temp\683472m3\nested\p1-g07d01\cellRegistered_fin.mat');
cregstruct_intermediates{3} = cell_registered_struct;

% meta reference map
meta_refmap = cregstruct_meta.cell_to_index_map;

% im footprints corrected
im_fp_corrected = cregstruct_meta.spatial_footprints_corrected;

% iterate through intermediate structures pulling footprints and indices
origin_fp_corrected = cell(1,20);
cell_to_index_map = zeros(size(meta_refmap,1), 20);


origin_sesh_ct = 0;
for im = 1:length(cregstruct_intermediates)

    % relevant cell nums
    all_im_row_nums = meta_refmap(meta_refmap(:,im)>0,im);

    if ~isempty(cregstruct_intermediates{im})
    
        for origin_sesh = 1:length(cregstruct_intermediates{im}.spatial_footprints_corrected)
            origin_sesh_ct = origin_sesh_ct+1;

            % fps
            origin_fp_corrected{origin_sesh_ct} = im_fp_corrected{im}(cregstruct_intermediates{im}.cell_to_index_map(:,origin_sesh)>0, :, :);
            
            % cell reg
            %
            
            % first add all cells from meta cell reg
            cell_to_index_map(:,origin_sesh_ct) = meta_refmap(:,im);
            
            % then replace cell numbers with rows in IM cell reg
            current_cell_idx = cell_to_index_map(:,origin_sesh_ct)>0;
            current_cell_nums = cell_to_index_map(current_cell_idx,origin_sesh_ct);
            cell_to_index_map(current_cell_idx,origin_sesh_ct) = cregstruct_intermediates{im}.cell_to_index_map(current_cell_nums,origin_sesh);

        end
    
    else
        origin_sesh_ct = origin_sesh_ct+1;
        
        origin_fp_corrected{origin_sesh_ct} = im_fp_corrected{im};
        cell_to_index_map(:,origin_sesh_ct) = meta_refmap(:,im);
        
    end
    
end


%
save_file_paths = ...
   [{'C:\Users\ampm1\Desktop\cellreg_temp\683472m3\spatial_footprints_nested_fin\spatial_footprints_probe02.mat'}...
    {'C:\Users\ampm1\Desktop\cellreg_temp\683472m3\spatial_footprints_nested_fin\spatial_footprints_probe03.mat'}...
    {'C:\Users\ampm1\Desktop\cellreg_temp\683472m3\spatial_footprints_nested_fin\spatial_footprints_probe04.mat'}...
    {'C:\Users\ampm1\Desktop\cellreg_temp\683472m3\spatial_footprints_nested_fin\spatial_footprints_problem_first_02.mat'}...
    {'C:\Users\ampm1\Desktop\cellreg_temp\683472m3\spatial_footprints_nested_fin\spatial_footprints_problem_first_03.mat'}...
    {'C:\Users\ampm1\Desktop\cellreg_temp\683472m3\spatial_footprints_nested_fin\spatial_footprints_problem_last_01.mat'}...
    {'C:\Users\ampm1\Desktop\cellreg_temp\683472m3\spatial_footprints_nested_fin\spatial_footprints_problem_last_02.mat'}...
    {'C:\Users\ampm1\Desktop\cellreg_temp\683472m3\spatial_footprints_nested_fin\spatial_footprints_problem_last_03.mat'}...
    {'C:\Users\ampm1\Desktop\cellreg_temp\683472m3\spatial_footprints_nested_fin\spatial_footprints_probe05.mat'}...
    {'C:\Users\ampm1\Desktop\cellreg_temp\683472m3\spatial_footprints_nested_fin\spatial_footprints_probe06.mat'}...
    {'C:\Users\ampm1\Desktop\cellreg_temp\683472m3\spatial_footprints_nested_fin\spatial_footprints_probe07.mat'}...
    {'C:\Users\ampm1\Desktop\cellreg_temp\683472m3\spatial_footprints_nested_fin\spatial_footprints_probe08.mat'}...
    {'C:\Users\ampm1\Desktop\cellreg_temp\683472m3\spatial_footprints_nested_fin\spatial_footprints_problem_first_05.mat'}...
    {'C:\Users\ampm1\Desktop\cellreg_temp\683472m3\spatial_footprints_nested_fin\spatial_footprints_problem_first_06.mat'}...
    {'C:\Users\ampm1\Desktop\cellreg_temp\683472m3\spatial_footprints_nested_fin\spatial_footprints_problem_last_04.mat'}...
    {'C:\Users\ampm1\Desktop\cellreg_temp\683472m3\spatial_footprints_nested_fin\spatial_footprints_problem_last_05.mat'}...
    {'C:\Users\ampm1\Desktop\cellreg_temp\683472m3\spatial_footprints_nested_fin\spatial_footprints_problem_last_06.mat'}...
    {'C:\Users\ampm1\Desktop\cellreg_temp\683472m3\spatial_footprints_nested_fin\spatial_footprints_probe01.mat'}...
    {'C:\Users\ampm1\Desktop\cellreg_temp\683472m3\spatial_footprints_nested_fin\spatial_footprints_problem_first_01.mat'}...
    {'C:\Users\ampm1\Desktop\cellreg_temp\683472m3\spatial_footprints_nested_fin\spatial_footprints_problem_first_04.mat'}];
%}
[~, sort_idx] = sort(save_file_paths);

%{
for isesh = 1:length(origin_fp_corrected)
    footprints = origin_fp_corrected{isesh};
    save(save_file_paths{isesh}, 'footprints');
end
%}

% sort
origin_fp_corrected = origin_fp_corrected(sort_idx);
cell_to_index_map = cell_to_index_map(:,sort_idx);

% build structure
clear cell_registered_struct
cell_registered_struct.spatial_footprints_corrected = origin_fp_corrected;
cell_registered_struct.cell_to_index_map = cell_to_index_map;



end


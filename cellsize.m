function mean_cell_sizes = cellsize(folderpath)

% get subject folders
subject_folders = get_folder_paths_all(folderpath, 0);

% preallcoate
mean_cell_sizes = nan(size(subject_folders,1),2);

% iterate through subjects
for isubj = 1:size(subject_folders,1)
isubj
    
    % get footprints folder
    footprint_folder_path = get_folder_paths_all(subject_folders{isubj}, 0);
    footprint_folder_path = footprint_folder_path(contains(footprint_folder_path, 'spatial_footprints'));

    % get footprint files
    footprint_files = get_file_paths_targeted(footprint_folder_path{1}, 'spatial_footprints');

    % iterate through sessions
    session_cellsizes = [];
    for isesh = 1:size(footprint_files,1)
        
        load(footprint_files{isesh}, 'footprints')

        for icell = 1:size(footprints,1)

            cell_ymax = max(sum(squeeze(footprints(icell,:,:)>0)));
            cell_xmax = max(sum(squeeze(footprints(icell,:,:)>0)'));

            session_cellsizes = [session_cellsizes; [cell_ymax cell_xmax]];

        end

    end
    
    % means
    cs_means = nanmean(session_cellsizes(session_cellsizes(:,1)>0 & session_cellsizes(:,2)>0,2));
    cs_medians = nanmedian(session_cellsizes(session_cellsizes(:,1)>0 & session_cellsizes(:,2)>0,2));
    
    
    
    mean_cell_sizes(isubj,1) = mean(cs_means);
    mean_cell_sizes(isubj,2) = mean(cs_medians);
    
end

end
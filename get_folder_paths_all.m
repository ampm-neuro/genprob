function item_paths = get_folder_paths_all(folderpath, depth, depth_count)
%iterates through each folder in path and returns path incl folder

%for only folders immediately in path, set depth to 0, do not provide an
%input for depth_count

item_paths = cell(1);

% control how many folders deep we go
if exist('depth_count', 'var')
    depth_count = depth_count+1;
    if exist('depth', 'var') && ~isnan(depth)
        if depth<depth_count
            return
        end
    else
        depth = nan;
    end

else
    depth_count = 0;
    
    if ~exist('depth', 'var')
        depth = nan;
    end
end


%items in path
file_list = dir(folderpath);
for iitem = 1:length(file_list)
    current_sesh = file_list(iitem).name;
    
    %omited folder and file names (folders named 'old' are invisible)
    if strcmp(current_sesh, '.') || strcmp(current_sesh, '..') || strcmp(current_sesh, 'old')
        continue
    end

    if isfolder([folderpath '\' file_list(iitem).name])
        item_paths = [item_paths; {[folderpath '\' file_list(iitem).name]}];
        item_paths_hold = get_folder_paths_all([folderpath '\' file_list(iitem).name], depth, depth_count);
        item_paths = [item_paths; item_paths_hold];
    end

end

%remove empty cells
item_paths = item_paths(find(~cellfun(@isempty, item_paths)));
end
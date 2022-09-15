function delete_empty_folders(folderpath)
% delete any empty folders

% get all folders
item_paths = get_folder_paths_all(folderpath);

% sort by length
[~,stringLength] = sort(cellfun(@length,item_paths),'descend');
item_paths = item_paths(stringLength);

% iterate through trying to delete (matlab doesnt allow deletion of folders
% unless they are empty
for i = 1:size(item_paths,1)
   [~] = rmdir(item_paths{i});
end

end


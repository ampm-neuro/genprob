function delete_todays_files(folderpath)
% deletes any files added today

%'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data'

% get all file paths
file_paths = get_file_paths_all(folderpath);
file_paths = file_paths(~contains(file_paths, 'old'));

% delete today's files
for iepath = 1:size(file_paths,1)
    
    % get file info
    file = dir(file_paths{iepath});
    mod_date = file.date(1:strfind(file.date, ' ')-1);
    
    % if modified date is today, delete
    if datetime(mod_date) == datetime(date)
        delete(file_paths{iepath})
    end
end


% get all folder paths
folder_paths = get_folder_paths_all(folderpath);
folder_paths = folder_paths(~contains(folder_paths, 'old'));

% delete today's files
for iepath = 1:size(folder_paths,1)
    
    % get file info
    fileinfo = getfileinfo(folder_paths{iepath});
    c_date = fileinfo.CreationDate;
    
    % if modified date is today, delete
    if datetime(c_date) == datetime(date)
        rmdir(folder_paths{iepath})
    end
end




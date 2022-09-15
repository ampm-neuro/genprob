function remove_subject(subject_id)
%moves all data files from a single subject into a folder (in cage file)
%called 'old'. Deletes empty session folders.

    %data folder
    umbrella_folder_path = 'E:\Projects\InProgress\GenProb\data\Gen_richards\repeat_probes\train_scrambled';

    %cage number
    cage_number = subject_id(1:end-2);

    % folder with cage
    umbrella_folder_path = [umbrella_folder_path '\' cage_number];

    %create 'old' storage folder if it doesnt already exist
    if ~exist([umbrella_folder_path '\old'], 'dir')
        mkdir([umbrella_folder_path '\old'])
    end

    %get relevant folder contents
    fld_cnts = folder_contents(umbrella_folder_path);

    %isolate folders
    folder_names_stages = {fld_cnts([fld_cnts(:).isdir]).name};

    %iterate through folders looking for subject_id files
    for ifile_stage = 1:length(folder_names_stages)

            %get relevant folder contents
            fld_cnts = folder_contents([umbrella_folder_path '\' folder_names_stages{ifile_stage}]);
            folder_names_sesh_nums = {fld_cnts([fld_cnts(:).isdir]).name};

            %iterate through folders looking for subject_id files
            for ifile_sesh = 1:length(folder_names_sesh_nums)
                
                %get relevant folder contents
                fld_cnts = folder_contents([umbrella_folder_path '\' folder_names_stages{ifile_stage} '\' folder_names_sesh_nums{ifile_sesh}]);
                file_names = {fld_cnts(~[fld_cnts(:).isdir]).name};

                %if subject file exists
                if isempty(file_names)

                    %delete folder
                    rmdir([umbrella_folder_path '\' folder_names_stages{ifile_stage} '\' folder_names_sesh_nums{ifile_sesh}]);

                elseif all(contains(file_names, subject_id))
                    
                    %move file
                    movefile([umbrella_folder_path '\' folder_names_stages{ifile_stage} '\' folder_names_sesh_nums{ifile_sesh} '\' file_names{contains(file_names, subject_id)}], [umbrella_folder_path '\old']);
                    %delete folder
                    rmdir([umbrella_folder_path '\' folder_names_stages{ifile_stage} '\' folder_names_sesh_nums{ifile_sesh}]);

                elseif any(contains(file_names, subject_id))

                    %move file
                    movefile([umbrella_folder_path '\' folder_names_stages{ifile_stage} '\' folder_names_sesh_nums{ifile_sesh} '\' file_names{contains(file_names, subject_id)}], [umbrella_folder_path '\old']);
                end
            end       
    end

end

function contents = folder_contents(fpath)

    %identify cage folder contents
    contents = dir(fpath);

    %delete shit folder things
    delete_idx = [];
    for iitem = 1:size(contents,1)
        if strcmp(contents(iitem).name, '.') || ...
                strcmp(contents(iitem).name, '..') || ...
                strcmp(contents(iitem).name, 'old')
            delete_idx = [delete_idx; iitem];
        end
    end
    contents(delete_idx) = [];

end


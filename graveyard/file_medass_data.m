function file_medass_data(folder_path)
%moves data files from an umbrella folder into subject>stage>session
%folders

%identify folder contents
immediate_files_and_folders = dir(folder_path);

%delete shit folder things
delete_idx = [];
for iitem = 1:size(immediate_files_and_folders,1)
    if strcmp(immediate_files_and_folders(iitem).name, '.') || ...
            strcmp(immediate_files_and_folders(iitem).name, '..') || ...
            strcmp(immediate_files_and_folders(iitem).name, 'old')
        delete_idx = [delete_idx; iitem];
    end
end
immediate_files_and_folders(delete_idx) = [];

%seperate files and folders
folder_names = {immediate_files_and_folders([immediate_files_and_folders(:).isdir]).name};
file_names = {immediate_files_and_folders(~[immediate_files_and_folders(:).isdir]).name};


for ifile = 1:length(file_names)
    
    %get relevant file info
    %
        %read output (table)
        fid = fopen([folder_path '\' file_names{ifile}]);
        medass_output = textscan(fid, '%s%s%s%s%s%s');

        %isolate session ID information
        end_id = find(strcmp(medass_output{1,1}, 'MSN:'));
        medass_info = [];
        for mi = 1:size(medass_output,2)
            medass_info = [medass_info medass_output{1,mi}(1:end_id)];
        end

        %mpc (medass program) file
        mpc_file = medass_info{10,2};
        str_start_flg = 'ampm_'; str_start = strfind(mpc_file, str_start_flg);
        str_end_flg = '_B'; str_end = strfind(mpc_file, str_end_flg);

        if ~isempty(str_end)
            mpc_file = mpc_file((str_start + length(str_start_flg)) : str_end - 1);
        else
            mpc_file = mpc_file((str_start + length(str_start_flg)) : end);
        end
        
        
        %EXCEPTIONS
        if strcmp(mpc_file, 'gen05_train02')
           mpc_file = 'gen06_train02'; 
        end

        %cage number
        cage_number = medass_info{4,2}(1:end-2);

        %subject (mouse) number
        mouse_num = medass_info{4,2}((end-1) : end);
        
    
    % build relevant folders
    %
        % ensure cage folder exists 
        if ~exist([folder_path '\' cage_number], 'dir')
            mkdir([folder_path '\' cage_number])
        end

        % ensure mpc-file folder exists
        if ~exist([folder_path '\' cage_number '\', mpc_file], 'dir')
            mkdir([folder_path '\' cage_number '\', mpc_file])
        end
        
        %check how many times the subject has run this mpc file
        all_paths = get_file_paths_all([folder_path '\' cage_number '\', mpc_file]);
        day = sum(contains(all_paths, [cage_number mouse_num])) + 1;
        day = num2str(day);
        if length(day)<2; day = ['0' day]; end
        
        % ensure day folder exists
        if ~exist([folder_path '\' cage_number '\', mpc_file, '\' day], 'dir')
            mkdir([folder_path '\' cage_number '\', mpc_file, '\' day])
        end
        
        
    % move files
        % move file
        fclose('all');
        movefile([folder_path '\' file_names{ifile}], [folder_path '\' cage_number '\', mpc_file, '\' day]);
    
    
end


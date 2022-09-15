function medass_loadandsave(folderpath)
% create mat files for each medassfile. then you can manually add them to
% behavioral data folder


% medass file paths
mfps = get_file_paths_all(folderpath);
mfps = mfps(~contains(mfps, '.mat'));

% load each file
for ipath = 1:size(mfps,1)
    
    % subject id
    subj_id = find_subj_id(mfps{ipath});
    
    % save mat file in same folder
    load_medass_save(mfps{ipath}, folderpath, subj_id);
    
end

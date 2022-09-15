function [moved_or_no] = check_m_ids(folderpath)


% known subject ids
subj_ids = {'644784m1','644784m2','644784m3','644784m4','644785m1',...
    '644785m2','644785m3','644785m4','644786m1','644786m2','644786m3'...
    '644786m4','645214m1','645214m2','645214m3','645214m4','645216m1'...
    '645216m2','645216m3','645216m4','647383m1','647383m2','647383m4'...
    '647718m1','647718m2','648721m1','648721m2','648721m3','648721m4'...
    '636838m1','636838m2'};

% get all file paths
file_paths = get_file_paths_all(folderpath);

% find paths that don't contain a known subject id
error_paths = file_paths(~contains(file_paths, subj_ids));

% move to error folder
for iepath = 1:size(error_paths,1)
    movefile(error_paths{iepath},'C:\Users\ampm1\Desktop\sandy_attach_error')
end


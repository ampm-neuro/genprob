function all_out = ALL(datafolder)
%generic function that runs through sessions

%should be input
folderpath = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\';
folderpath = [folderpath datafolder];

%preallocate
%
all_out = [];


%subject folders
file_list_subjects = dir(folderpath);
file_list_subjects(1:2) = [];
for isubject = 1:length(file_list_subjects)
    current_subj = file_list_subjects(isubject).name;
    
    %session files
    file_list_sessions = dir([folderpath '\' current_subj]);
    file_list_sessions(1:2) = [];
    for isession = 1 : length(file_list_sessions)
        current_sesh = file_list_sessions(isession).name;

        %load data file
        load([folderpath '\' current_subj '\' current_sesh])

        %load output
        all_out = [all_out; ];

        
    end
end


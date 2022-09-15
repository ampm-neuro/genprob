function [nld_out, all_subjects] = num_learning_days(datafolder, training_stages)
%compute learning days for each subject on each stage


% learning folderpath
folderpath = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\';
folderpath = [folderpath datafolder];

% unique subjects (cell)
all_subjects = struct2cell(dir(folderpath));
all_subjects = all_subjects(1,3:end)';
all_subjects(contains(all_subjects, 'old')) = [];

% preallocate output
nld_out = nan(size(all_subjects,1), length(training_stages));

% iterate through subjects
for isubj = 1:size(all_subjects,1)

    % subject path
    folderpath_subj = [folderpath all_subjects{isubj}];
    
    % iterate through training stages
    istage_ct = 0;
    for istage = training_stages
        istage_ct = istage_ct+1;
        
        stage_str = num2str(istage);
        if length(stage_str)<2
            stage_str = ['0' stage_str];
        end
        
        % all learning sessions
        relevant_session_files = path_comp(folderpath_subj, ['train' stage_str]);
        num_session_files = length(relevant_session_files);
        
        % load output
        nld_out(isubj, istage_ct) = num_session_files;
    
    end
end
function [session_files_out, included_subjects] = get_file_paths_targeted(string_patterns, day_delay, new_learn_num)
%returns list of session files from specific experimental conditions that
%match the input strings in the cell 'string_pattern'


%input checks
if size(day_delay,1)>size(day_delay,2)
    day_delay = day_delay';
end
if size(new_learn_num,1)>size(new_learn_num,2)
    new_learn_num = new_learn_num';
end

%where files are stored
data_files_folder = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\Gen_richards';

%get files containing all string patterns
all_files = get_file_paths(data_files_folder)';
str_constrained_files = all_files;
for istr = 1:length(string_patterns)
    str_constrained_files = str_constrained_files(contains(str_constrained_files, string_patterns{istr}));
end

%unique subjects
unq_subj = [];
for ipf = 1:size(str_constrained_files,2)
    start_pt = length(data_files_folder)+2;  
    end_pt = strfind(str_constrained_files{ipf}(start_pt:end), '\'); 
        end_pt = start_pt + end_pt(1) - 2;
    subj_text = str_constrained_files{ipf}(start_pt:end_pt);
    unq_subj{ipf} = subj_text;
end
unq_subj = unique(unq_subj)';


%only include files from subjects that experienced the input delay day condition(s)
if ~isempty(day_delay)
    delay_num_subjs = [];
    for iunq = 1:length(unq_subj)
        for idd = day_delay

            day_delay_str = num2str(idd); 
            if length(day_delay_str)<2
                day_delay_str = ['0' day_delay_str];
            end

            subj_idx = contains(all_files, unq_subj{iunq});
            delay_num_idx = contains(all_files, [day_delay_str 'd']);

            if ~isempty(all_files(subj_idx & delay_num_idx))
                delay_num_subjs = [delay_num_subjs; unq_subj(iunq)];
            end

        end
    end
    %reduce
    if isempty(delay_num_subjs)
        disp('requested sessions do not exist');
        session_files_out=[];
        included_subjects=[];
        return
    end
    str_constrained_files = str_constrained_files(contains(str_constrained_files, delay_num_subjs));
end

%only include files from subjects that experienced the input new learning condition(s)
if ~isempty(new_learn_num)
    NL_num_subjs = [];
    for iunq = 1:length(unq_subj)
        for idd = new_learn_num

            new_learn_str = num2str(idd); 
            if length(new_learn_str)<2
                new_learn_str = ['0' new_learn_str]; 
            end

            subj_idx = contains(all_files, unq_subj{iunq});
            NL_num_idx = contains(all_files, ['newlearn' new_learn_str]);

            if ~isempty(all_files(subj_idx & NL_num_idx))
                NL_num_subjs = [NL_num_subjs; unq_subj(iunq)];
            end

        end
    end

    %reduce
    if isempty(NL_num_subjs)
        display('requested sessions do not exist');
        session_files_out=[];
        included_subjects=[];
        return
    end
    str_constrained_files = str_constrained_files(contains(str_constrained_files, NL_num_subjs));
    if isempty(str_constrained_files)
        display('requested sessions do not exist');
        session_files_out=[];
        included_subjects=[];
        return
    end
end


% find remaining subjects
all_subjects = [];
for ipf = 1:size(str_constrained_files,2)
    start_pt = length(data_files_folder)+2;  
    end_pt = strfind(str_constrained_files{ipf}(start_pt:end), '\'); 
        end_pt = start_pt + end_pt(1) - 2;
    subj_text = str_constrained_files{ipf}(start_pt:end_pt);
    all_subjects{ipf} = subj_text;
end

%output
session_files_out = str_constrained_files';
included_subjects = unique(all_subjects)';









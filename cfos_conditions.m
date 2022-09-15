% identify condition for each subject in a cfos folder

% conditions
% 1 = mevar, probe 02
% 2 = mevar, probe 07
% 3 = mevar, probe 08
% 4 = hivar, probe 02
% 5 = hivar, probe 07
% 6 = hivar, probe 08

load('cfos_cts.mat')

exp_folder_path = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone';
exp_folders = get_folder_paths_all(exp_folder_path, 0);
cfos_folders = exp_folders(contains(exp_folders, 'cfos'));

subject_ids = [];
subject_condition = [];
subj_coefs = [];
for icond = 1:size(cfos_folders,1)
    
    if icond == 1
        icond_modifier = 3;
    else
        icond_modifier = 0;
    end
    
    subject_folders = get_folder_paths_all(cfos_folders{icond}, 0);
    
    for isubj = 1:size(subject_folders,1)
        
        % subject ids
        last_slash = strfind(subject_folders{isubj}, '\'); last_slash = last_slash(end);
        subj_id = subject_folders{isubj}(last_slash+1:end);
        subject_ids = [subject_ids; subj_id];
        
        
        % subject conditions
        probe_session_files = [get_file_paths_targeted(subject_folders{isubj}, 'preprobe'); ...
            get_file_paths_targeted(subject_folders{isubj}, 'postprobe')];
               
        if contains(probe_session_files{end}, 'postprobe_02')
            subject_condition = [subject_condition; 3+icond_modifier];
        elseif contains(probe_session_files{end}, 'postprobe_01')
            subject_condition = [subject_condition; 2+icond_modifier];
        elseif contains(probe_session_files{end}, 'preprobe_02')
            subject_condition = [subject_condition; 1+icond_modifier];
        else
            error('weird')
        end
        
        
        % coefs of last probe
        load(probe_session_files{end}, 'trl_mtx')
        wait_times = wait_times_prep(trl_mtx, 1, 4);
        [~, coefEsts, modelFun] = ampm_normal_logistic_fit_algo(1:41, wait_times, [nanmean(wait_times) 20 0 1]);
        %[~, coefEsts, modelFun] = ampm_normal_logistic_fit(1:41, wait_times, [nanmean(wait_times) 20 0 1]);
        subj_coefs = [subj_coefs; coefEsts];
        
    end
    

end

[subject_ids, sort_idx] = sortrows(subject_ids);
subject_condition = subject_condition(sort_idx);
subj_coefs = subj_coefs(sort_idx,:);






function [wait_durations, unq_frq, p_dist, day_bin_fielddiff_zs] = GenLearn(datafolder, sesh)
%plot wait time curves over all learning days

% folderpath
folderpath = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\';
folderpath = [folderpath datafolder];


%subject folders
file_list_subjects = dir(folderpath);
file_list_subjects(1:2) = [];
all_trl_mtx = [];
day_bin_fielddiff_zs = [];

for isubject = 1:length(file_list_subjects)
    current_subj = file_list_subjects(isubject).name;

    %session files
    file_list_sessions = dir([folderpath '\' current_subj]);
    file_list_sessions(1:2) = [];


    for isession = sesh

        %skip if session file does not exist
        if length(file_list_sessions)<isession
            continue
        end

        current_sesh = file_list_sessions(isession).name;

        %load data file
        load([folderpath '\' current_subj '\' current_sesh])

        
            %set trial minimum
            %{
            min_trials = 50;
            if size(trl_mtx,1)<min_trials
                continue
            end 
            %}
            %limit trials
            %{
            max_trials = 100;
            if size(trl_mtx,1)>max_trials
                trl_mtx = trl_mtx(1:max_trials,:);
            end
            %}
            
            
         
        if isession == 9
            
            %load z difference
            if size(trl_mtx,1)>50
                [~, wait_durations_holdz] = wait_times(trl_mtx(1:50,:), medass_cell, 0);
            else
                [~, wait_durations_holdz] = wait_times(trl_mtx, medass_cell, 0);
            end
        else
                    
            %load z difference
            [~, wait_durations_holdz] = wait_times(trl_mtx, medass_cell, 0);
            
        end

            in_field_waits = cell2mat(wait_durations_holdz(p_dist>.5));
            out_field_waits = cell2mat(wait_durations_holdz(p_dist<.5));

            field_wait_diff = nanmean(in_field_waits) - nanmean(out_field_waits);
            pool_std = mean([nanstd(in_field_waits) nanstd(out_field_waits)]);
            day_bin_fielddiff_zs = [day_bin_fielddiff_zs; [isubject field_wait_diff/pool_std]];
        
        %load output
        all_trl_mtx = [all_trl_mtx; trl_mtx];

    end
    
end

%wait times for current session
[~, wait_durations, ~, unq_frq, ~] = wait_times(all_trl_mtx, medass_cell, 0);











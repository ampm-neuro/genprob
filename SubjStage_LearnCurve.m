function [rich_cell, poor_cell, rwd_recency_bias_cell] = SubjStage_LearnCurve(folder_path, subject, stage)
%computes wait times for rich and poor sessions within a single training
%stage for a single subject


%subject folderpath
subj_folder_path = [folder_path '\' subject];

%all files
item_paths = get_file_paths_all(subj_folder_path);

%only keep files from the desired stage
if ~isstring(stage); stage = num2str(stage); end
if length(stage)==1; stage = ['0' stage]; end
item_paths = item_paths(contains(item_paths, ['novar' stage]) | contains(item_paths, ['lovar' stage]) | contains(item_paths, ['mevar' stage]) | contains(item_paths, ['hivar' stage])  | contains(item_paths, ['exvar' stage]) | contains(item_paths, ['ctl' stage])  | contains(item_paths, ['exvar' stage]));


%preallocate output cells
rich_cell = cell(1, length(item_paths));
poor_cell = cell(1, length(item_paths));
rwd_recency_bias_cell = cell(1, length(item_paths));

%iterate though paths
for isession = 1:length(item_paths)
    
    %load session
    load(item_paths{isession})
    
    %compute waits
    [wait_durations, wd_freq] = wait_times_prep(trl_mtx,2);
    [prob_dist, pd_freq] = rwd_prob_by_freq(medass_cell);
    
    %load cells
    rich_cell{isession} = wait_durations(ismember(wd_freq,pd_freq(prob_dist>0.5)));
    poor_cell{isession} = wait_durations(ismember(wd_freq,pd_freq(prob_dist<0.5)));
    rwd_recency_bias_cell{isession} = rwd_recency_dep(trl_mtx, 1);
    
end


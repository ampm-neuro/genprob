function [subject_cell, rwd_bias_subj_cell] = ALL_stage_LearnCurves(folderpath, stage)
%plots learning curves from each training stage for all animals. input
%dictates how to divide up the sessions. e.g., if num_stage_samples == 3,
%then the total number of sessions will be divided into early, middle, and
%late, and plotted as 3 data points per subject.

%subject folders
file_list_subjects = dir(folderpath);
file_list_subjects(1:2) = [];

%iterate through subjects
subj_count = 0;
for isubject = 1:length(file_list_subjects)
    current_subj = file_list_subjects(isubject).name;
    subj_count = subj_count +1;

    %compute wait times over stage for this subject
    [rich_cell, poor_cell, rwd_recency_bias_cell] = SubjStage_LearnCurve(folderpath, current_subj, stage);
    subject_cell{subj_count} = [rich_cell; poor_cell];
    rwd_bias_subj_cell{subj_count} = rwd_recency_bias_cell;
    
end



function ALL_opto_discrim_full(problems, folder_details)
% plot tone discrimination with and without opto


%% find unique subjects

% base folder path
fp = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\';
%fp = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_optoExp_acc\';

train_mevar_optoExp_' region '\'];
fp2 = ['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_hivar_optoExp_' region '\'];

session_files = get_file_paths_targeted(fp1, {'mevar0', 'opto-'});
session_files = [session_files; get_file_paths_targeted(fp2, {'hivar0', 'opto-'})];
% find unique subjects
subjs = [];
for isf = 1:size(session_files,1)
    subjs = [subjs; {session_files{isf}(strfind(session_files{isf}, '\\')+2 : strfind(session_files{isf}, '\\')+9)}];
end
subjs = unique(subjs);


%% preallocate
ALL_discrims_optoOFF = 

%% iterate through each problem




%% stats
%{
% pairwise
[~,p_first,~,stats_first] = ttest(prefs_for_rich_means{1}, prefs_for_rich_means{2});
[~,p_last,~,stats_last] = ttest(prefs_for_rich_means{3}, prefs_for_rich_means{4});

% mixed model
tbl = table(prefs_for_rich_means(:), all_means_sesh_id(:), all_means_opto_id(:), all_means_subj_id(:), all_means_problem_id(:),...
    'VariableNames',{'pref','sesh','opto','subj','problem'});
tbl.sesh = categorical(tbl.sesh);
tbl.opto = categorical(tbl.opto);
tbl.subject = categorical(tbl.subj);
tbl.problem = categorical(tbl.problem);

model_str = 'pref~opto+sesh+opto*sesh+(1|subject)+(1|problem)';
lme = fitlme(tbl, model_str);
%}










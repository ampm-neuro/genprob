function probe_vs_behavior_opto
% correlate measures of probe behavior with problem behavior

%% load data
% bar charts of problem performance
[prob_all_means, prob_prefs_for_rich, prob_rich_waits, prob_poor_waits, prob_sesh_id,...
    prob_opto_id, prob_problem_id, prob_subj_ids] = ALL_opto_discrim_full(1:6); % Input desired session numbers
prob_subj_ids = cell2mat(prob_subj_ids);

% probe over probe change
[probe_brON, probe_bpON, probe_brOFF, probe_bpOFF, probe_subj_ids, probe_problemNums] = probe_change_opto(1:7);

%% common subjects only

% common subjects
[shared_subjs] = intersect(prob_subj_ids,probe_subj_ids, 'rows', 'stable');
subjidx_problem = ismember(prob_subj_ids, shared_subjs, 'rows');
subjidx_probe = ismember(probe_subj_ids, shared_subjs, 'rows');

    % problem data
    %prob_all_means = prob_all_means(subjidx_problem,:);
    prob_prefs_for_rich = prob_prefs_for_rich(subjidx_problem,:);
    prob_rich_waits = prob_rich_waits(subjidx_problem,:);
    prob_poor_waits = prob_poor_waits(subjidx_problem,:);
    prob_sesh_id = prob_sesh_id(subjidx_problem,:);
    prob_opto_id = prob_opto_id(subjidx_problem,:);
    prob_problem_id = prob_problem_id(subjidx_problem,:);
    prob_subj_ids = prob_subj_ids(subjidx_problem,:);
    
    % probe data
    probe_brON = probe_brON(subjidx_probe,:);
    probe_bpON = probe_bpON(subjidx_probe,:);
    probe_brOFF = probe_brOFF(subjidx_probe,:);
    probe_bpOFF = probe_bpOFF(subjidx_probe,:);
    probe_subj_ids = probe_subj_ids(subjidx_probe,:);
    probe_problemNums = probe_problemNums(subjidx_probe,:);
    
    
%% sort match subjects/session pairs

% problem pairs
[prob_subj_ids_srt, prob_subj_ids_srt_idx] = sortrows(prob_subj_ids); 
[~,~,prob_subject_id_nums] = unique(prob_subj_ids_srt, 'rows', 'stable');
unsrt_index(prob_subj_ids_srt_idx) = 1:length(prob_subject_id_nums);
prob_subject_id_nums = prob_subject_id_nums(unsrt_index);

% probe pairs
[probe_subj_ids_srt, probe_subj_ids_srt_idx] = sortrows(probe_subj_ids); 
[~,~,probe_subject_id_nums] = unique(probe_subj_ids_srt, 'rows', 'stable');
unsrt_index(probe_subj_ids_srt_idx) = 1:length(probe_subject_id_nums);
probe_subject_id_nums = probe_subject_id_nums(unsrt_index);

% matches
prob_subjProb = [prob_subject_id_nums prob_problem_id(:,1)];
probe_subjProb = [probe_subject_id_nums probe_problemNums];
[shared_subjsProbs] = intersect(prob_subjProb,probe_subjProb, 'rows', 'stable');
subjProbidx_problem = ismember(prob_subjProb, shared_subjsProbs, 'rows');
subjProbidx_probe = ismember(probe_subjProb, shared_subjsProbs, 'rows');

    % problem data
    %prob_all_means = prob_all_means(subjidx_problem,:);
    prob_prefs_for_rich = prob_prefs_for_rich(subjProbidx_problem,:);
    prob_rich_waits = prob_rich_waits(subjProbidx_problem,:);
    prob_poor_waits = prob_poor_waits(subjProbidx_problem,:);
    prob_sesh_id = prob_sesh_id(subjProbidx_problem,:);
    prob_opto_id = prob_opto_id(subjProbidx_problem,:);
    prob_problem_id = prob_problem_id(subjProbidx_problem,:);
    prob_subj_ids = prob_subj_ids(subjProbidx_problem,:);
    prob_subjProb = prob_subjProb(subjProbidx_problem,:);
    
    % probe data
    probe_brON = probe_brON(subjProbidx_probe,:);
    probe_bpON = probe_bpON(subjProbidx_probe,:);
    probe_brOFF = probe_brOFF(subjProbidx_probe,:);
    probe_bpOFF = probe_bpOFF(subjProbidx_probe,:);
    probe_subj_ids = probe_subj_ids(subjProbidx_probe,:);
    probe_problemNums = probe_problemNums(subjProbidx_probe,:);
    probe_subjProb = probe_subjProb(subjProbidx_probe,:);


%% plot 

%
% drop outlier
laserAffect_richtonechange = probe_brON-probe_brOFF;
del_idx = laserAffect_richtonechange>.15;
% problem data
prob_prefs_for_rich(del_idx,:) = [];
prob_rich_waits(del_idx,:) = [];
prob_poor_waits(del_idx,:) = [];
prob_sesh_id(del_idx,:) = [];
prob_opto_id(del_idx,:) = [];
prob_problem_id(del_idx,:) = [];
prob_subj_ids(del_idx,:) = [];
prob_subjProb(del_idx,:) = [];
% probe data
probe_brON(del_idx,:) = [];
probe_bpON(del_idx,:) = [];
probe_brOFF(del_idx,:) = [];
probe_bpOFF(del_idx,:) = [];
probe_subj_ids(del_idx,:) = [];
probe_problemNums(del_idx,:) = [];
probe_subjProb(del_idx,:) = [];
%}

% WAIT TIME
% if laser similarly affects RICH wait problem behavior and probe
laserAffect_richwaits = prob_rich_waits(:,4)-prob_rich_waits(:,3);    
laserAffect_richtonechange = probe_brON-probe_brOFF;
figure;
plot(laserAffect_richwaits, laserAffect_richtonechange, 'o', 'color', .7.*[1 1 1])
[r, p] = fit_line(laserAffect_richwaits, laserAffect_richtonechange)
xlabel('Laser affect on Rich Tone Wait Time')
ylabel('Laser affect on Rich Probe change')
title(['r=' num2str(r) '; p=' num2str(p)])

% if laser similarly affects POOR wait problem behavior and probe
laserAffect_poorwaits = prob_poor_waits(:,2)-prob_poor_waits(:,1);
laserAffect_poortonechange = probe_bpON-probe_bpOFF;
figure;
plot(laserAffect_poorwaits, laserAffect_poortonechange, 'o', 'color', .7.*[1 1 1])
[r, p] = fit_line(laserAffect_poorwaits, laserAffect_poortonechange)
xlabel('Laser affect on Poor Tone Wait Time')
ylabel('Laser affect on Poor Probe change')
title(['r=' num2str(r) '; p=' num2str(p)])



% PREFERENCE
% if laser similarly affects RICH preference problem behavior and probe
laserAffect_richpref = prob_prefs_for_rich(:,4)-prob_prefs_for_rich(:,3);
laserAffect_richtonechange = probe_brON-probe_brOFF;
figure;
plot(laserAffect_richpref, laserAffect_richtonechange, 'o', 'color', .7.*[1 1 1])
[r, p] = fit_line(laserAffect_richpref, laserAffect_richtonechange)
xlabel('Laser affect on Rich Tone Pref')
ylabel('Laser affect on Rich Probe change')
title(['r=' num2str(r) '; p=' num2str(p)])

% if laser similarly affects POOR preference problem behavior and probe
laserAffect_poorpref = prob_prefs_for_rich(:,2)-prob_prefs_for_rich(:,1);
laserAffect_poortonechange = probe_bpON-probe_bpOFF;
figure;
plot(laserAffect_poorpref, laserAffect_poortonechange, 'o', 'color', .7.*[1 1 1])
[r, p] = fit_line(laserAffect_poorpref, laserAffect_poortonechange)
xlabel('Laser affect on Poor Tone Pref')
ylabel('Laser affect on Poor Probe change')
title(['r=' num2str(r) '; p=' num2str(p)])



% CROSS (poor affects vs rich affects)
figure;
plot(laserAffect_poorwaits, laserAffect_richwaits, 'o', 'color', .7.*[1 1 1])
[r, p] = fit_line(laserAffect_poorwaits, laserAffect_richwaits)
xlabel('Laser affect on Poor Tone Wait Time')
ylabel('Laser affect on Rich Tone Wait Time')
title(['r=' num2str(r) '; p=' num2str(p)])

figure;
plot(laserAffect_poorpref, laserAffect_richpref, 'o', 'color', .7.*[1 1 1])
[r, p] = fit_line(laserAffect_poorpref, laserAffect_richpref)
xlabel('Laser affect on Poor Tone Pref')
ylabel('Laser affect on Rich Tone Pref')
title(['r=' num2str(r) '; p=' num2str(p)])

figure;
plot(laserAffect_poortonechange, laserAffect_richtonechange, 'o', 'color', .7.*[1 1 1])
[r, p] = fit_line(laserAffect_poortonechange, laserAffect_richtonechange)
xlabel('Laser affect on Poor Probe change')
ylabel('Laser affect on Rich Probe change')
title(['r=' num2str(r) '; p=' num2str(p)])


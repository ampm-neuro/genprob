function [weights, feds, num_trials, performances] = ALL_weight
% plots correlation between number of trials completed and weight as a
% percentage of FF

%should be input
folderpath = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\Gen';

%preallocate
weights = [];
feds = [];
num_trials = [];
performances = [];

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

        %performance stats
        num_trials = [num_trials; size(trl_mtx,1) ];
        %performances = [performances; ];
        
        if exist('weight', 'var')
            weights = [weights; [weight isubject isession]];
            feds = [feds; [fed isubject isession]];
        else
            weights = [weights; [nan isubject isession]];
            feds = [feds; [nan isubject isession]];
        end
        
        clearvars weight fed
        
    end
end


%figures

%weight (today)
figure; hold on
[r, p] = fit_line(weights(:,1).*0.8, num_trials);
ylabel('number of trials')
xlabel('weight (% of FF)')
title([num2str(r) '; ' num2str(p)])
for i = unique(weights(:,2))'
   plot(weights(weights(:,2)==i, 1).*0.8,  num_trials(weights(:,2)==i), 'o')
end

%feed (yesterday)
figure; hold on
feed_nm1 = feds(:, 1);
[r, p] = fit_line(feed_nm1, num_trials);
ylabel('number of trials')
xlabel('yesterday feed (g)')
title([num2str(r) '; ' num2str(p)])
for i = unique(weights(:,2))'
   plot(feed_nm1(feds(:,2)==i,1),  num_trials(feds(:,2)==i), 'o')
end

%feed (two days ago)
figure; hold on
ylabel('number of trials')
xlabel('2days ago feed (g)')
feed_hold = [];
numtrl_hold = [];
for i = unique(weights(:,2))'
    feed_nm2 = feds(feds(:,2)==i,1);
    feed_nm2 = feed_nm2(1:end-1, 1);
    numtrials_nm2 = num_trials(feds(:,2)==i);
    numtrials_nm2 = numtrials_nm2(2:end);
    plot(feed_nm2, numtrials_nm2, 'o')
    feed_hold = [feed_hold; feed_nm2];
    numtrl_hold = [numtrl_hold; numtrials_nm2];
end
[r, p] = fit_line(feed_hold, numtrl_hold);
title([num2str(r) '; ' num2str(p)])
for i = unique(weights(:,2))'
    feed_nm2 = feds(feds(:,2)==i,1);
    feed_nm2 = feed_nm2(1:end-1, 1);
    numtrials_nm2 = num_trials(feds(:,2)==i);
    numtrials_nm2 = numtrials_nm2(2:end);
    plot(feed_nm2, numtrials_nm2, 'o')
    feed_hold = [feed_hold; feed_nm2];
    numtrl_hold = [numtrl_hold; numtrials_nm2];
end

%feed (3 days ago)
n_var = 2;
figure; hold on
ylabel('number of trials')
xlabel([num2str(n_var) 'days ago feed (g)'])
feed_hold = [];
numtrl_hold = [];
for i = unique(weights(:,2))'
    feed_nm2 = feds(feds(:,2)==i,1);
    feed_nm2 = feed_nm2(1:end-n_var, 1);
    numtrials_nm2 = num_trials(feds(:,2)==i);
    numtrials_nm2 = numtrials_nm2(n_var+1:end);
    plot(feed_nm2, numtrials_nm2, 'o')
    feed_hold = [feed_hold; feed_nm2];
    numtrl_hold = [numtrl_hold; numtrials_nm2];
end
[r, p] = fit_line(feed_hold, numtrl_hold);
title([num2str(r) '; ' num2str(p)])
for i = unique(weights(:,2))'
    feed_nm2 = feds(feds(:,2)==i,1);
    feed_nm2 = feed_nm2(1:end-n_var, 1);
    numtrials_nm2 = num_trials(feds(:,2)==i);
    numtrials_nm2 = numtrials_nm2(n_var+1:end);
    plot(feed_nm2, numtrials_nm2, 'o')
    feed_hold = [feed_hold; feed_nm2];
    numtrl_hold = [numtrl_hold; numtrials_nm2];
end


%session number
figure; hold on
[r, p] = fit_line(feds(:,3), num_trials);
ylabel('number of trials')
xlabel('session number')
title([num2str(r) '; ' num2str(p)])
for i = unique(weights(:,2))'
   plot(feds(feds(:,2)==i, 3),  num_trials(feds(:,2)==i), 'o')
end


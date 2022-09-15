function all_out = ALL_inoutfield(sessions)
%generic function that runs through sessions

%should be input
folderpath = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\Gen_richards\';

%preallocate
%
all_out = [];


%subject folders
file_list_subjects = dir(folderpath);
file_list_subjects(1:2) = [];
for isubject = 1:1:length(file_list_subjects)
    current_subj = file_list_subjects(isubject).name;
    
    %session files
    file_list_sessions = dir([folderpath '\' current_subj]);
    file_list_sessions(1:2) = [];
    
    for isession = sessions
        
        %skip if session file does not exist
        if length(file_list_sessions)<isession
            continue
        end
        
        current_sesh = file_list_sessions(isession).name;

        %load data file
        load([folderpath '\' current_subj '\' current_sesh])

        %compute wait durations
        [mean_wait_durations, ~, ~, unq_frq, p_dist] = wait_times(trl_mtx,medass_cell,0);

        %load output
        all_out = [all_out; [nanmean(mean_wait_durations(p_dist>.5)) nanmean(mean_wait_durations(p_dist<.5))]];

        
    end
end

%plot
celle{1} = all_out(:,1); celle{2} = all_out(:,2);
errorbar_plot(celle, 1)
ylim([0 25])

%paired ttest
[~, p, ~, stats] = ttest(all_out(:,1), all_out(:,2));

%add stats to plot
hold on;
plot([1 2], [1 1].*(max(all_out(:))*1.15), 'k-')
if p < 0.05
    text(1.5, max(all_out(:))*1.25, '*', 'color', [0 0 0], 'fontsize', 25)
else
    text(1.42, max(all_out(:))*1.27, 'n.s.', 'color', [0 0 0], 'fontsize', 17)
end
title(['t(' num2str(stats.df) ')=' num2str(stats.tstat) ', p=' num2str(p)])

%beautification
ylabel('Wait times (s)')
xticklabels({'In field', 'Out of field'})
xticks([1 2])
    




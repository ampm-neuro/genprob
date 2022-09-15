function [subj_pass, testpass, tstats, pvals, subjects_sessions] = ALL_testpass(datafolder, varargin)
%generic function that runs through sessions


%default is each subjects last session, otherwise input desired session
if nargin == 2
    sessions_hold = varargin{1};
else
    sessions_hold = [];
end

% folderpath
folderpath = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\';
folderpath = [folderpath datafolder];

%preallocate
%
tstats = [];
pvals = [];
testpass = [];
subjects = [];
current_sessions= [];

%subject folders
file_list_subjects = dir(folderpath);
file_list_subjects(1:2) = [];
for isubject = 1:1:length(file_list_subjects)
    current_subj = file_list_subjects(isubject).name;
    
    if contains(current_subj, 'old')
        continue
    end
    
    %session files
    file_list_sessions = dir([folderpath '\' current_subj]);
    file_list_sessions(1:2) = [];
    
    if ~isempty(sessions_hold)
        sessions = sessions_hold;
    else
        sessions = length(file_list_sessions);
    end
    
    for isession = sessions
        
        %skip if session file does not exist
        if length(file_list_sessions)<isession
            continue
        end
        
        current_sesh = file_list_sessions(isession).name;
        
        % exceptions for probe sessions
        if contains(current_sesh, 'probe') || contains(current_sesh, 'remind')
            tstats = [tstats; [nan nan]];
            pvals = [pvals; [nan nan nan]];
            testpass = [testpass; true];
            subjects = [subjects; current_subj];
            current_sessions = [current_sessions; 'probe_session!'];
            continue
            
        elseif contains(current_sesh, 'new')
            
        end

        %load data file
        load([folderpath '\' current_subj '\' current_sesh])

        %compute wait durations
        [wait_durations_all, wda_freq] = wait_times_prep(trl_mtx,2);
        unq_frqs = unique(wda_freq);

        %test differences between rich and poor tone wait times
        
            %anova
            anova_grp = []; 
            for itone = 1:length(unq_frqs)
                anova_grp = [ anova_grp; itone*ones( size( wait_durations_all(wda_freq==unq_frqs(itone)) ) ) ];
            end 
            anova_p = anovan(wait_durations_all, anova_grp, 'display', 'off');

            %posthoc ttests
            [~, pval1, ~, stats1] = ttest2(wait_durations_all(wda_freq==unq_frqs(2)), wait_durations_all(wda_freq==unq_frqs(1)));
            [~, pval2, ~, stats2] = ttest2(wait_durations_all(wda_freq==unq_frqs(3)), wait_durations_all(wda_freq==unq_frqs(1)));
            tstat = [stats1.tstat stats2.tstat];
            

        %load output
        tstats = [tstats; tstat];
        pvals = [pvals; [pval1 pval2 anova_p]];
        testpass = [testpass; tstat(1)<0 & tstat(2)<0 & pval1<0.05 & pval2<0.05 & anova_p<0.05 & length(wait_durations_all(wda_freq==unq_frqs(1)))>=2];
        subjects = [subjects; current_subj];
        current_sessions = [current_sessions; current_sesh(end-13:end)];
        
        %exception for New Learning sessions    
        if contains(current_sesh, 'new')
            if contains(current_sesh, '14')
                testpass(end) = true;
            else
                testpass(end) = false;
            end
        end
        
    end
    
    sessions = [];
end

%summary output
subj_pass = subjects(logical(testpass),:);
subjects_sessions = [subjects repmat(' : ', size(subjects,1),1) current_sessions];



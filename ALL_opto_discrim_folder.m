function ALL_opto_discrim_folder(problems, folder_detail)
% plot tone discrimination with and without opto for each problem input for
% all subjects in a given folder


%% find unique subjects

% base folder path
fp = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\';
fp = [fp folder_detail '\'];

% all opto session files in folder
session_files = get_file_paths_targeted(fp, {'opto-'})

% find unique subjects among session files
subjects = [];
for isf = 1:size(session_files,1)
    subjects = [subjects; {session_files{isf}(strfind(session_files{isf}, '\\')+2 : strfind(session_files{isf}, '\\')+9)}];
end
subjects = unique(subjects);



%% iterate through each problem
problem_cell = cell(2,length(problems));
prob_ct = 0;
for iproblem = problems
    prob_ct = prob_ct+1;
    problem_cell{prob_ct} = cell(1,2);
    
    
    % iterate first and last
    for first_last = 1:2
    
        % compute and plot discriminations on this problem
        [discrims_optoOFF, discrims_optoON, session_files] = opto_discrim_multisubj(subjects, iproblem, first_last, fp);
    
        % load
         problem_cell{prob_ct}{first_last}{1} = discrims_optoOFF;
         problem_cell{prob_ct}{first_last}{2} = discrims_optoON;
         problem_cell{prob_ct}{first_last}{3} = [];
             
            for isesh = 1:size(session_files,1)
                 
                 subj = find_subj_id(session_files{isesh});
                 if isempty(subj)
                    subj = 'no file ';
                 end
                 problem_cell{prob_ct}{first_last}{3} = [problem_cell{prob_ct}{first_last}{3}; subj];

            end
             
            %problem_cell{prob_ct}{first_last}{1}
            %problem_cell{prob_ct}{first_last}{2}
            %problem_cell{prob_ct}{first_last}{3}
          
    end
    
end


%% match subjects

% unique subjects
all_subjects = [];
for iproblem = 1:length(problems)
    for first_last = 1:2
        all_subjects = [all_subjects;  problem_cell{iproblem}{first_last}{3}];
    end
end
unique_subjs = unique(all_subjects, 'rows');
unique_subjs = unique_subjs(~ismember(unique_subjs, 'no file ', 'rows'),:);
unique_subjs = sortrows(unique_subjs);

% align sessions with unique subject list
for iproblem = 1:length(problems)
    for first_last = 1:2
        
        % align subject list
        [~,subj_idx] = ismember(problem_cell{iproblem}{first_last}{3}, unique_subjs, 'rows');
        subj_idx(subj_idx==0) = [];
        problem_cell{iproblem}{first_last}{3} = repmat('no file ', size(unique_subjs,1),1);
        problem_cell{iproblem}{first_last}{3}(subj_idx,:) = unique_subjs(subj_idx,:);
        
        % align data
        hold_first = problem_cell{iproblem}{first_last}{1};
            problem_cell{iproblem}{first_last}{1} = nan(size(unique_subjs,1),1);
            problem_cell{iproblem}{first_last}{1}(subj_idx,:) = hold_first(subj_idx);
        hold_second = problem_cell{iproblem}{first_last}{2};
            problem_cell{iproblem}{first_last}{2} = nan(size(unique_subjs,1),1);
            problem_cell{iproblem}{first_last}{2}(subj_idx,:) = hold_second(subj_idx);
        
    end
end



%% plots

% BAR
for iproblem = 1:length(problems)
    figure; hold on
    
    colors = distinguishable_colors(size(unique_subjs,1));
    legend_trick(colors, '-')
    
    xhold = [1 2; 3.5 4.5];
    for first_last = 1:2
        errorbar_barplot(problem_cell{iproblem}{first_last}(1:2), 1, xhold(first_last,:), colors);
    end
    
    % aesthetics
    xlim([0 6])
    ylim([-4 6])
    xticks(sort(xhold(:)))
    xticklabels({'FirstOFF', 'FirstON', 'LastOFF', 'LastON'})
    legend(unique_subjs, 'location', 'northeastoutside')
    ylabel('Preference for rich tone (CohenD)')
    
end


% LINE
xhold = [(.8:1:5.8)' (1.2:1:6.2)'];
y_color = [0.9882 0.7294 0.0118];
b_color = [0.05 0.29 0.65];
figure; hold on
for iproblem = 1:length(problems)
    
    %colors = distinguishable_colors(size(unique_subjs,1));
    %legend_trick(colors, '-')
    
    
    plot_line_OFF = [problem_cell{iproblem}{1}(1) problem_cell{iproblem}{2}(1)];
    plot_line_ON = [problem_cell{iproblem}{1}(2) problem_cell{iproblem}{2}(2)];
    errorbar_plot(plot_line_OFF, 1, xhold(iproblem,:), repmat(.8.*y_color, size(unique_subjs,1),1), y_color);
    errorbar_plot(plot_line_ON, 1, xhold(iproblem,:), repmat(.8.*b_color, size(unique_subjs,1),1), b_color);
    
    %
    %for first_last = 1:2
    %    errorbar_plot(problem_cell{iproblem}{first_last}(1:2), 1, xhold(first_last,:), colors);
    %end
    
    % aesthetics
    xlim([0 ceil(max(xhold(iproblem,:)))])
    ylim([-4 6])
    xticks(sort(mean(xhold,2)))
    xticklabels({'First problem', 'Last problem'})
    %legend(unique_subjs, 'location', 'northeastoutside')
    ylabel('Preference for rich tone (CohenD)')
    plot(xlim, [1 1].*0, 'k--')

    
end





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










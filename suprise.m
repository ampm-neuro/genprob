function [suprise_cell, discrim_cell, suprise_mtx_interp, problem_paths] = suprise(subject_ids)

% plot rich-tone preference overlaid with suprise
% subject_ids = unique_subjs('\train_mevar_imaging_hpc\');


%% basic info
wdw_size = 10;
load('unqfrq41', 'unqfrq41')

% behavioral folder
beh_fld = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\';


% get all first-day paths for each subject
problem_paths = cell(1, size(subject_ids,1));
for isubj = 1:size(subject_ids,1)
    all_subj_paths = get_file_paths_targeted(beh_fld, {subject_ids(isubj,:)});
    for iproblem = 1:6
        all_problem_files = [all_subj_paths(contains(all_subj_paths, ['var0' num2str(iproblem)])); ...
            all_subj_paths(contains(all_subj_paths, ['ctl0' num2str(iproblem)]))];
        problem_paths{isubj} = [problem_paths{isubj}; all_problem_files(1)];
    end
end
  
    
%% iterate through learning sessions getting rich-tone preferences
all_subj_wait_prefs = cell(1, size(subject_ids,1));
all_subj_problem_idx = cell(1, size(subject_ids,1));
for isubj = 1:length(problem_paths)
    
    % iterate through first-day sessions
    for isesh = 1:size(problem_paths{isubj},1)
        % compute trial-by-trial preference for rich tone
        [wait_prefs_sesh] = continuous_learning_singlesesh_suprise(problem_paths{isubj}{isesh}, wdw_size);
        all_subj_wait_prefs{isubj}{isesh} = wait_prefs_sesh;
        all_subj_problem_idx{isubj}{isesh} = repmat(isesh, size(wait_prefs_sesh));
    end
    
end

discrim_cell = all_subj_wait_prefs;



%% iterate through learning sessions identifying suprise
all_subj_suprise = cell(1, size(subject_ids,1));
all_subj_suprise_inst = cell(1, size(subject_ids,1));
for isubj = 1:length(problem_paths)
    
    % iterate through first-day sessions
    for isesh = 1:size(problem_paths{isubj},1)
        
        % iterate through trials
        load(problem_paths{isubj}{isesh}, 'trl_mtx')
        all_subj_suprise{isubj}{isesh} = nan(size(trl_mtx,1),1);
        all_subj_suprise_inst{isubj}{isesh} = nan(size(trl_mtx,1),1);
        
        % rich tone
        rich_tone = mode(floor(trl_mtx(trl_mtx(:,3)==1,2)));
        rich_idx = floor(trl_mtx(:,2))==rich_tone;
        %reward idx
        rwd_idx = ~isnan(trl_mtx(:,11));
        
        % suprise idx
        all_subj_suprise{isubj}{isesh}(rich_idx & rwd_idx) = 1; % yes rich, yes rwd
        all_subj_suprise{isubj}{isesh}(~rich_idx & ~rwd_idx) = 2; % no rich, no rwd
        all_subj_suprise{isubj}{isesh}(rich_idx & ~rwd_idx) = 3; % yes rich, no rwd
        all_subj_suprise{isubj}{isesh}(~rich_idx & rwd_idx) = 4; % no rich, yes rwd
        
        % compute proportion of suprising trials in window
        for itrl = 1 : size(trl_mtx,1)
            local_window = (itrl-wdw_size : itrl+wdw_size); % all past trials 
            %local_window = (1 : itrl); % all past trials 
            %local_window = (itrl+1:size(trl_mtx,1)); % future trials (not sig)
            local_window = local_window(local_window>0 & local_window<=size(trl_mtx,1));
            
            unpredicted_events = sum(ismember(all_subj_suprise{isubj}{isesh}(local_window), [3 4]));
            predicted_events = sum(ismember(all_subj_suprise{isubj}{isesh}(local_window), [1 2]));
            
            all_subj_suprise_inst{isubj}{isesh}(itrl) = (unpredicted_events-predicted_events)/(unpredicted_events+predicted_events);
        end
        
        % multiply suprise by subjects instantaneous prediction
        %all_subj_suprise_inst_fin{isubj}{isesh} = all_subj_suprise_inst{isubj}{isesh};
        all_subj_suprise_inst_fin{isubj}{isesh} = all_subj_suprise_inst{isubj}{isesh}.*all_subj_wait_prefs{isubj}{isesh};
                %all_subj_suprise_inst_fin{isubj}{isesh} = all_subj_suprise_inst{isubj}{isesh}.*(all_subj_wait_prefs{isubj}{isesh}./abs(all_subj_wait_prefs{isubj}{isesh}));

    end
    
end

suprise_cell = all_subj_suprise_inst_fin;



%% subject plots
%{

for isubj = 1:length(problem_paths)

    % subject plot
    figure; hold on;
    
    % plot preferences
    preference_vect = cell2mat(all_subj_wait_prefs{isubj}');
    plot(1:length(preference_vect), preference_vect, '-', 'color', 0.5.*[1 1 1])

    % plot suprise
    suprise_vect = cell2mat(all_subj_suprise_inst_fin{isubj}');
    plot(1:length(suprise_vect), suprise_vect, 'r-')
        
    % problem boundaries
    ylim([-1 1]); xlim([1 length(suprise_vect)])
    problem_vect = cell2mat(all_subj_problem_idx{isubj}');
    for isesh = unique(problem_vect)'
        plot((find(problem_vect==isesh, 1, 'last')+0.5).*[1 1], ylim, 'k--')
    end
    
    % no preference
    plot(xlim, [0 0], 'k--');
    drawnow
end
%}


%% overall plots
%
all_subj_wait_prefs_interp = all_subj_wait_prefs;
all_subj_suprise_inst_interp = all_subj_suprise_inst_fin;
all_pref_mtx = [];
all_supr_mtx = [];
for isubj = 1:length(problem_paths)
    for isesh = 1:length(problem_paths{isubj})
   
        % warp preferences
        preference_vect = all_subj_wait_prefs_interp{isubj}{isesh};
        all_subj_wait_prefs_interp{isubj}{isesh} = interp1(1:length(preference_vect), preference_vect, linspace(1,length(preference_vect),100));

        % warp suprise
        suprise_vect = all_subj_suprise_inst_interp{isubj}{isesh};
        all_subj_suprise_inst_interp{isubj}{isesh} = interp1(1:length(suprise_vect), suprise_vect, linspace(1,length(suprise_vect),100));
        
    end
    
    % cell2mat
    all_pref_mtx = [all_pref_mtx; cell2mat(all_subj_wait_prefs_interp{isubj})];
    all_supr_mtx = [all_supr_mtx; cell2mat(all_subj_suprise_inst_interp{isubj})];

end

% out
suprise_mtx_interp = all_supr_mtx;

% means and stds
all_pref_mean = nanmean(all_pref_mtx);
all_pref_ste = nanstd(all_pref_mtx)./sqrt(sum(~isnan(all_pref_mtx)));
all_supr_mean = nanmean(all_supr_mtx);
all_supr_ste = nanstd(all_supr_mtx)./sqrt(sum(~isnan(all_supr_mtx)));
%
figure; hold on
for isubj = 1:size(all_pref_mtx,1)
    
    plot(all_supr_mtx(isubj,:), 'color', [1 .85 .85], 'linewidth', 1)
  
end
%
    plot(all_pref_mean, 'color', .4.*[1 1 1], 'linewidth', 2)
        plot(all_pref_mean-all_pref_ste, 'color', .65.*[1 1 1], 'linewidth', 1)
        plot(all_pref_mean+all_pref_ste, 'color', .65.*[1 1 1], 'linewidth', 1)  
    plot(all_supr_mean, 'color', [0.7 0 0], 'linewidth', 2)
        plot(all_supr_mean-all_supr_ste, 'color', [1 0 0], 'linewidth', 1)
        plot(all_supr_mean+all_supr_ste, 'color', [1 0 0], 'linewidth', 1)




% problem bounds
ylim([-1 1]); set(gca,'TickLength',[0, 0]);
plot([100.5:100:600]'.*[1 1], ylim, 'k--')
% no preference
plot(xlim, [0 0], 'k--');
%}













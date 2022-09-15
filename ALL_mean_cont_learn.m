function [all_prefs, subj_idx, problem_idx, sesh_idx] = ALL_mean_cont_learn(training_type)
% single line for all subjects showing mean behavioral curve

% window size
wdw_size = 11;

% all behavioral files
training_group = 'train_mevar_meta';
session_paths01 = get_file_paths_targeted(....
    ['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\' training_type '\'], 'var0');
session_paths02 = get_file_paths_targeted(....
    ['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\' training_type '\'], 'ctl');
session_paths = sort([session_paths01; session_paths02]);

% preallocate
all_prefs = [];
subj_idx = [];
problem_idx = [];
prog_idx = [];


% find unique subjects
unq_subjects = [];
for isf = 1:size(session_paths,1)
    local_path = session_paths{isf};
    search_str = '\\';
    rng_vals = [2 9];
    unq_subjects = [unq_subjects; {local_path(strfind(local_path, search_str)+rng_vals(1) : strfind(local_path, search_str)+rng_vals(2))}];
end
unq_subjects = unique(unq_subjects);

% iterate through subjects
for isubj = 1:size(unq_subjects,1)
    
    % sessions paths featuring subject
    subj_session_paths = session_paths(contains(session_paths, unq_subjects{isubj}));

    
    % unique training problems
    unq_trainprobs = [];
    for isf = 1:size(subj_session_paths,1)
        local_path = subj_session_paths{isf};
        search_str = unq_subjects{isubj};
        rng_vals = [length(search_str)+1 length(search_str)+5];
        unq_trainprobs = [unq_trainprobs; {local_path(strfind(local_path, search_str)+rng_vals(1) : strfind(local_path, search_str)+rng_vals(2))}];
    end
    unq_trainprobs = unique(unq_trainprobs);
    
    

    % iterate through training problems (typically 1:6)
    for itp = 1:size(unq_trainprobs,1)
        
        
        % preallocate
        prob_prefs = [];
        
        % subject session paths featuring training problem
        trainprob_session_paths = subj_session_paths(contains(subj_session_paths, unq_trainprobs{itp}));
        
        % iterate through each session
        pass_fail = 0;
        for isesh = 1:size(trainprob_session_paths,1)
            if pass_fail == 1; continue; end

            % compute tone preferences
            [tone_prefs] = continuous_learning_singlesesh(trainprob_session_paths{isesh}, wdw_size);
            
            % interp to 100 trials
            old_x = (1:length(tone_prefs))./length(tone_prefs);
            new_x = linspace(min(old_x), max(old_x), 100);
            tone_prefs = interp1(old_x, tone_prefs, new_x)';
            
            % report if cont learning function produces nans
            if isnan(sum(tone_prefs))
                trainprob_session_paths{isesh}
            end
                
            
            % load session preferences
            prob_prefs = [prob_prefs; tone_prefs];
            
            % was this the crit session?
            load(trainprob_session_paths{isesh})
            [wait_durations_all, wda_freq] = wait_times_prep(trl_mtx,2);
            [prob_dist, pd_freq] = rwd_prob_by_freq(medass_cell);
            rich_freq = pd_freq(prob_dist>0.5);
            poor_freqs = pd_freq(prob_dist<0.5);
            [~, pval, ~, stats] = ttest2(wait_durations_all(wda_freq==poor_freqs), wait_durations_all(wda_freq==rich_freq));
            tstat = stats.tstat;
            if tstat<0 && pval<=0.01
                pass_fail = 1;
            else
                pass_fail = 0;
            end

        end
        
        % load overall
        all_prefs = [all_prefs; prob_prefs];
        prog_idx = [prog_idx; ((1:length(prob_prefs))./length(prob_prefs))'];
        subj_idx = [subj_idx; repmat(isubj, size(prob_prefs))];
        problem_idx = [problem_idx; repmat(itp, size(prob_prefs))];
        
    end
end

% binning correction
prog_idx(prog_idx==0) = realmin;

% bining
slide_window_size = 0.025;
HiLo_bins = linspace(0,1,(1/slide_window_size) + 1);
lobins = HiLo_bins(1:end-1);
hibins = HiLo_bins(2:end);

% figure
h1 = figure; hold on
for iprob = 1:6
    
    local_progs = prog_idx(problem_idx==iprob);
    local_prefs = all_prefs(problem_idx==iprob);
    
    % compute means and stes for each bin
    prob_means = [];
    prob_stds = [];
    for ibin = 1:length(lobins)
        
        prob_means = [prob_means; nanmean(local_prefs(local_progs>lobins(ibin) & local_progs<=hibins(ibin)))];
        prob_stds = [prob_stds; nanstd(local_prefs(local_progs>lobins(ibin) & local_progs<=hibins(ibin)))];

    end
    
    % plot
    figure(h1)
    plot(((1:length(prob_means))./length(prob_means))+iprob, prob_means, 'k-', 'linewidth', 2)
    plot(((1:length(prob_means))./length(prob_means))+iprob, prob_means+prob_stds, 'k-', 'linewidth', 1)
    plot(((1:length(prob_means))./length(prob_means))+iprob, prob_means-prob_stds, 'k-', 'linewidth', 1)
    plot(iprob.*[1 1], [-1 1], 'k-')
    
    
end
    
% aesthetics
hold on; plot(xlim, [0 0], 'k--')
set(gca,'TickLength',[0, 0]); box off;
ylim([-.4 .8])
ylabel('Preference for rich tone')
xlabel('Training problem')
    
    
    
    
    
    
    
    
    
    
function [all_prefs, subj_idx, problem_idx, first_mtx, last_mtx, all_sesh_ct, learn_rate_mtx, delta_mtx, min_discrim_mtx, unq_subjects, first_mtx_rich, first_mtx_poor, last_mtx_rich, last_mtx_poor] = ALL_cont_learn_delta_dot(plot_on, training_type)
% single line for all subjects showing mean behavioral curve

% plot on?
if nargin==0
    plot_on = 0;
end

% training_group = 'train_hivar';
    

% window size
wdw_size = 11;

% CHOSE PROPORTION OF LEARNING CURVE HERE
%
% delta proportion of curve (<=0.5) 
delta_prop = 0.20;

% all behavioral files
%fldr = 'train_mevar_meta';
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
all_smooth_rich_waits = [];
all_smooth_poor_waits = [];

% find unique subjects
unq_subjects = [];
for isf = 1:size(session_paths,1)
    local_path = session_paths{isf};
    search_str = '\\';
    rng_vals = [2 9];
    unq_subjects = [unq_subjects; {local_path(strfind(local_path, search_str)+rng_vals(1) : strfind(local_path, search_str)+rng_vals(2))}];
end
unq_subjects = unique(unq_subjects);

% hardish preallocate
all_sesh_ct = nan(length(unq_subjects),6);

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
        prob_rich_waits = [];
        prob_poor_waits = [];
        
        % subject session paths featuring training problem
        trainprob_session_paths = subj_session_paths(contains(subj_session_paths, unq_trainprobs{itp}));
        
        % iterate through each session
        pass_fail = 0;
        for isesh = 1:size(trainprob_session_paths,1)
            if pass_fail == 1; continue; end

            % compute tone preferences
            [tone_prefs, rich_waits, poor_waits] = continuous_learning_singlesesh(trainprob_session_paths{isesh}, wdw_size);
            
            % interp to 100 trials
            %{
            old_x = (1:length(tone_prefs))./length(tone_prefs);
            new_x = linspace(min(old_x), max(old_x), 100);
            tone_prefs = interp1(old_x, tone_prefs, new_x)';
            %}
            
            % report if cont learning function produces nans
            if isnan(sum(tone_prefs))
                trainprob_session_paths{isesh}
            end
                
            
            % load session preferences
            prob_rich_waits = [prob_rich_waits; rich_waits];
            prob_poor_waits = [prob_poor_waits; poor_waits];
            prob_prefs = [prob_prefs; tone_prefs];
            if isnan(all_sesh_ct(isubj,itp))
                all_sesh_ct(isubj,itp) = length(tone_prefs);
            else
                all_sesh_ct(isubj,itp) = all_sesh_ct(isubj,itp)+length(tone_prefs); 
            end
            
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
        all_smooth_rich_waits = [all_smooth_rich_waits; prob_rich_waits];
        all_smooth_poor_waits = [all_smooth_poor_waits; prob_poor_waits];
        
    end
end

% compute days to crit and total learning
first_mtx = nan(size(all_sesh_ct));
last_mtx = nan(size(all_sesh_ct));
min_discrim_mtx = nan(size(all_sesh_ct));
for isubj = 1:size(unq_subjects,1)
    si = subj_idx==isubj;
    for iprob = 1:6
        pi = problem_idx==iprob;
                
        % subject and problem preferences
        local_progs = prog_idx(si & pi);
        local_prefs = all_prefs(si & pi);
        local_rich_waits = all_smooth_rich_waits(si & pi);
        local_poor_waits = all_smooth_poor_waits(si & pi);
        
        if ~isempty(local_prefs) && ~isempty(local_prefs)
        
            % compute preference
            first_prop = mean(local_prefs(local_progs<=delta_prop));
            first_rich = mean(local_rich_waits(local_progs<=delta_prop));
            first_poor = mean(local_poor_waits(local_progs<=delta_prop));
            last_prop = mean(local_prefs(local_progs>=(1-delta_prop)));
            last_rich = mean(local_rich_waits(local_progs>=(1-delta_prop)));
            last_poor = mean(local_poor_waits(local_progs>=(1-delta_prop)));

            % load
            first_mtx(isubj, iprob) = first_prop;
            first_mtx_rich(isubj, iprob) = first_rich;
            first_mtx_poor(isubj, iprob) = first_poor;
            last_mtx(isubj, iprob) = last_prop;
            last_mtx_rich(isubj, iprob) = last_rich;
            last_mtx_poor(isubj, iprob) = last_poor;

            if first_prop<-.4
               fm =  unq_subjects{isubj};
            end
            if first_prop<0 & last_prop<0
               bm =  unq_subjects{isubj};
            end
            
            % num trials to certain level of discrim (try .2)
            %{
            min_discrim = 0.3;
            try
                min_discrim_mtx(isubj, iprob) = min(local_progs(local_prefs>=min_discrim));
            catch
            end
            %}

        end
        
    end
end


% FIGURE ISSUE
first_mtx(first_mtx<-.4) = nan;


% ONLY COMPLETE MICE
%{
nan_idx = isnan(sum(first_mtx,2)) | isnan(sum(last_mtx,2));
first_mtx = first_mtx(~nan_idx,:);
last_mtx = last_mtx(~nan_idx,:);
all_sesh_ct = all_sesh_ct(~nan_idx,:);
min_discrim_mtx = min_discrim_mtx(~nan_idx,:);
unq_subjects = unq_subjects(~nan_idx);
first_mtx_rich = first_mtx_rich(~nan_idx,:);
first_mtx_poor = first_mtx_poor(~nan_idx,:);
last_mtx_rich = last_mtx_rich(~nan_idx,:);
last_mtx_poor = last_mtx_poor(~nan_idx,:);
%}

% delta matrix
delta_mtx = last_mtx-first_mtx;

% learning rate matrix
learn_rate_mtx = delta_mtx./all_sesh_ct;


% mixed model prep
subj_num_mtx = repmat((1:size(first_mtx,1))', 1, size(first_mtx,2)); % subject matrix
stage_num_mtx = repmat(1:size(first_mtx,2), size(first_mtx,1), 1); % stage matrix
model_str = 'dv~ProbeBeh+(1|subject)+(ProbeBeh-1|subject)+(1|problem)+(ProbeBeh-1|problem)';

if plot_on == 1
    % plot
    %colors = bone(15); colors = colors(7:13,:);
    colors = bone(10); colors = colors(2:7,:);

    figure; hold on; 
    for iprob = 1:6
        plot(first_mtx(:, iprob), all_sesh_ct(:, iprob), 'o', 'color', colors(iprob,:))
        plot(nanmean(first_mtx(:, iprob)), nanmean(all_sesh_ct(:, iprob)), '.', 'color', colors(iprob,:), 'markersize', 50)
    end
    [r p] = fit_line(first_mtx(:), all_sesh_ct(:));

            % mixed model
            tbl = table(all_sesh_ct(:), first_mtx(:), subj_num_mtx(:), stage_num_mtx(:),'VariableNames',{'dv','ProbeBeh','subject','problem'});
            tbl.subject = categorical(tbl.subject);
            tbl.problem = categorical(tbl.problem);
            lme_NumTrials = fitlme(tbl,model_str)

    % aesthetics
    set(gca,'TickLength',[0, 0]); box off;
    %ylim([0 900])
    hold on; plot([0 0], ylim, 'k--')
    xlim([-.3 .9])
    xlabel('Early discrimination')
    ylabel('Number of training trials required to reach criterion')
    title(['early discrim vs num trials; r=' num2str(r) ', p=' num2str(p)])

    figure; 
    for iprob=1:6

        errorbar_plot([{first_mtx(:, iprob)}, {last_mtx(:, iprob)}], 1, [iprob-.25 iprob+.25])

    end
    ylim([-.3 .91])
    yticks(-.3:.1:.9)
    xlim([0.5 iprob+0.5])
    hold on; plot(xlim, [0 0], 'k--')

    figure; hold on; 
    for iprob = 1:6
        plot(first_mtx(:, iprob), last_mtx(:, iprob), 'o', 'color', colors(iprob,:))
        plot(nanmean(first_mtx(:, iprob)), nanmean(last_mtx(:, iprob)), '.', 'color', colors(iprob,:), 'markersize', 50)
    end
    [r p] = fit_line(first_mtx(:), last_mtx(:));

            % mixed model
            tbl = table(last_mtx(:), first_mtx(:), subj_num_mtx(:), stage_num_mtx(:),'VariableNames',{'dv','ProbeBeh','subject','problem'});
            tbl.subject = categorical(tbl.subject);
            tbl.problem = categorical(tbl.problem);
            lme_LateDiscr = fitlme(tbl,model_str)

    set(gca,'TickLength',[0, 0]); box off;
    title(['early discrim vs late discrim; r=' num2str(r) ', p=' num2str(p)])


    %{
    figure; [r p] = fit_line(first_mtx(:), learn_rate_mtx(:))
    axis_hold = axis;
    figure; hold on; 
    for iprob = 1:6
        plot(first_mtx(:, iprob), learn_rate_mtx(:, iprob), 'o', 'color', colors(iprob,:))
        plot(nanmean(first_mtx(:, iprob)), nanmean(learn_rate_mtx(:, iprob)), '.', 'color', colors(iprob,:), 'markersize', 50)
    end
    axis(axis_hold);
    set(gca,'TickLength',[0, 0]); box off;
    title('First vs LearnRate')
    %}


    %first_mtx(isnan(sum(first_mtx,2)),:) = nan;
    %{
    figure; errorbar_plot([{first_mtx(:, 1)}, {first_mtx(:, 2)},...
        {first_mtx(:, 3)}, {first_mtx(:, 4)}, {first_mtx(:, 5)}, {first_mtx(:, 6)}], 1)
    title('Early discrimination')
    %}

    %{
    [~, twop] = ttest(first_mtx(:, 1), first_mtx(:, 2))
    [~, threep] = ttest(first_mtx(:, 1), first_mtx(:, 3))
    [~, fourp] = ttest(first_mtx(:, 1), first_mtx(:, 4))
    [~, fivep] = ttest(first_mtx(:, 1), first_mtx(:, 5))
    [~, sixp] = ttest(first_mtx(:, 1), first_mtx(:, 6))

    [~, one_to_0] = ttest(first_mtx(:, 1))
    [~, two_to_0] = ttest(first_mtx(:, 2))
    [~, three_to_0] = ttest(first_mtx(:, 3))
    %}
    %{
    figure; errorbar_plot([{learn_rate_mtx(:, 1)}, {learn_rate_mtx(:, 2)},...
        {learn_rate_mtx(:, 3)}, {learn_rate_mtx(:, 4)}, {learn_rate_mtx(:, 5)}, {learn_rate_mtx(:, 6)}])

    figure; errorbar_plot([{delta_mtx(:, 1)}, {delta_mtx(:, 2)},...
        {delta_mtx(:, 3)}, {delta_mtx(:, 4)}, {delta_mtx(:, 5)}, {delta_mtx(:, 6)}])


    figure; errorbar_plot([{min_discrim_mtx(:, 1)}, {min_discrim_mtx(:, 2)},...
        {min_discrim_mtx(:, 3)}, {min_discrim_mtx(:, 4)}, {min_discrim_mtx(:, 5)}, {min_discrim_mtx(:, 6)}])

    %}
end
    
    
    
    
    
    
    
    
    
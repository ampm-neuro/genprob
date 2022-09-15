function [poor_tone_rvals, rich_tone_rvals, poor_tone_waits, rich_tone_waits, rich_poor_tone_nums] = image_corr_popCorr_waits(subject, target_session, comparison_sessions, cell_regist_mtx, tses)
% breaks wait times from current probe into thirds (low tones,
% (mevar-rewarded) middle tones, and high tones) and correlates them with
% the average population distance from the target probe to the two input 
% comparison probes



%% compute each trial response for each neuron in each target session
targ_mtx = cell(1, length(target_session));
rich_poor_tone_nums = cell(1,2);
for itarg = 1:length(target_session)
    % load session files
    [trl_mtx_targ, trl_idx_targ, frame_times_targ, traces_targ] = load_trl_mtx(subject, target_session(itarg));
    
    % included trials
    incl_trials = 1:size(trl_mtx_targ,1);
    incl_trials = incl_trials(ismember(incl_trials,unique(trl_idx_targ)) & (trl_mtx_targ(:,3)==0)');
    
    % compute activity matrix
    [~, ~, ~, ~, ~, targ_mtx{itarg}] = image_mean_activity_timewarp_noRWD(trl_mtx_targ, trl_idx_targ, frame_times_targ, traces_targ, 1:size(traces_targ,1), incl_trials, tses);
    close;close
    
    % tone numbers
    load('unqfrq41', 'unqfrq41')
    all_tones = unique(floor(trl_mtx_targ(:,2)));
    all_tone_nums = find(ismember(unqfrq41,all_tones));
    rich_poor_tone_nums{1} = all_tone_nums(all_tone_nums>=16 & all_tone_nums<=26);
    rich_poor_tone_nums{2} = all_tone_nums(all_tone_nums<16 | all_tone_nums>26);
end



%% compute mean trial response for each neuron in each comparison session
comp_mtx = cell(1, length(comparison_sessions));
for icomp = 1:length(comparison_sessions)
    
    % load session files
    [trl_mtx, trl_idx, frame_times, traces] = load_trl_mtx(subject, comparison_sessions(icomp));

    % compute activity matrix
    [~, ~, ~, ~, comp_mtx{icomp}] = image_mean_activity_timewarp_noRWD(trl_mtx, trl_idx, frame_times, traces, 1:size(traces,1), unique(trl_idx), tses);
    close;close
end



%% compute wait times
probe_freqs_waits = trl_mtx_targ(incl_trials, [2 12]);
probe_freqs_waits(:,1) = floor(probe_freqs_waits(:,1));



%% compute average correlation
all_trial_correlations = nan(size(probe_freqs_waits,1), length(comparison_sessions));
for itrl = 1:size(probe_freqs_waits,1)
    
    % compute target session trial matrix
    trl_actmtx = nan(length(targ_mtx{1}), size(targ_mtx{1}{1},2));
    for ineuron = 1:length(targ_mtx{1})
        trl_actmtx(ineuron,:) = targ_mtx{1}{ineuron}(itrl,:);
    end

    
    % compare to each comp session
    common_neuron_rs = nan(length(targ_mtx{1}),length(comparison_sessions));
    for icomp = 1:length(comparison_sessions)
        
        % only shared neurons
        common_idx = cell_regist_mtx(:,target_session)>0 & cell_regist_mtx(:,comparison_sessions(icomp))>0;
        targ_trl_mtx_common = trl_actmtx(cell_regist_mtx(common_idx,target_session),:);
        comp_trl_mtx_common = comp_mtx{icomp}(cell_regist_mtx(common_idx,comparison_sessions(icomp)),:);
        
        % iterate through neurons
        for ineuron = 1:size(targ_trl_mtx_common,1)
            common_neuron_rs(ineuron, icomp) = corr(comp_trl_mtx_common(ineuron,:)', targ_trl_mtx_common(ineuron,:)');
        end
    end
    
    %load
    all_trial_correlations(itrl,:) = nanmean(common_neuron_rs);
    
end



%% find difference between comparisons
if length(comparison_sessions)==2
    all_trial_correlations = all_trial_correlations(:,2) - all_trial_correlations(:,1);
end



%% relevant tone values

% tone freqs
load('unqfrq41', 'unqfrq41')

% tone numbers thirds (inclusive)
lo_tone_nums = 1:15;
lo_tone_freqs = unqfrq41(lo_tone_nums);
med_tone_nums = 16:26;
med_tone_freqs = unqfrq41(med_tone_nums);
hi_tone_nums = 27:41;
hi_tone_freqs = unqfrq41(hi_tone_nums);



%% plot probe wait times for target and comparison sessions

% target
for itarg = target_session
    [trl_mtx] = load_trl_mtx(subject, itarg);
    figure; hold on
    wait_times_plot(trl_mtx,3, 0.8.*[1 1 1]);
    wait_times_plot_tonecolor(trl_mtx);
    title(['Session ' num2str(itarg)])
    ylim([0 60])
    plot(mean([lo_tone_freqs(end) med_tone_freqs(1)]).*[1 1], ylim, 'k--')
    plot(mean([med_tone_freqs(end) hi_tone_freqs(1)]).*[1 1], ylim, 'k--')
end

% comparison
for icomp = comparison_sessions
    [trl_mtx] = load_trl_mtx(subject, icomp);
    figure; hold on
    wait_times_plot(trl_mtx,3, 0.8.*[1 1 1]);
    wait_times_plot_tonecolor(trl_mtx);
    title(['Session ' num2str(icomp)])
    ylim([0 60])
    plot(mean([lo_tone_freqs(end) med_tone_freqs(1)]).*[1 1], ylim, 'k--')
    plot(mean([med_tone_freqs(end) hi_tone_freqs(1)]).*[1 1], ylim, 'k--')
end



%% correlate pop distances and wait times

% plot correlations
figure
tcolors = parula(71); tcolors = tcolors(26:66,:);

% low
subplot(1,3,1); hold on
trial_corrs_lo = all_trial_correlations(ismember(probe_freqs_waits(:,1), lo_tone_freqs));
wait_times_lo = probe_freqs_waits(ismember(probe_freqs_waits(:,1), lo_tone_freqs), 2);
if ~isempty(wait_times_lo)
    trialTones_lo = probe_freqs_waits(ismember(probe_freqs_waits(:,1), lo_tone_freqs), :);
    [r,p] = fit_line(trial_corrs_lo, wait_times_lo);
    for idot = lo_tone_nums
        tone_idx = trialTones_lo(:,1) == unqfrq41(idot);
        plot(trial_corrs_lo(tone_idx),  trialTones_lo(tone_idx, 2), 'o', 'color', tcolors(idot,:));
    end
    title(['r=' num2str(round(r*100)/100) ';p=' num2str(round(p*1000)/1000)])
end
    ylim([0 60]); set(gca,'TickLength',[0, 0]); box off;


% low
subplot(1,3,2); hold on
trial_corrs_med = all_trial_correlations(ismember(probe_freqs_waits(:,1), med_tone_freqs));
wait_times_med = probe_freqs_waits(ismember(probe_freqs_waits(:,1), med_tone_freqs), 2);
if ~isempty(wait_times_med)
    trialTones_med = probe_freqs_waits(ismember(probe_freqs_waits(:,1), med_tone_freqs), :);
    [r,p] = fit_line(trial_corrs_med, wait_times_med);
    for idot = med_tone_nums
        tone_idx = trialTones_med(:,1) == unqfrq41(idot);
        plot(trial_corrs_med(tone_idx),  trialTones_med(tone_idx, 2), 'o', 'color', tcolors(idot,:));
    end
    title(['r=' num2str(round(r*100)/100) ';p=' num2str(round(p*1000)/1000)])
end
    ylim([0 60]); set(gca,'TickLength',[0, 0]); box off;

% low
subplot(1,3,3); hold on
trial_corrs_hi = all_trial_correlations(ismember(probe_freqs_waits(:,1), hi_tone_freqs));
wait_times_hi = probe_freqs_waits(ismember(probe_freqs_waits(:,1), hi_tone_freqs), 2);
if ~isempty(wait_times_hi)
    trialTones_hi = probe_freqs_waits(ismember(probe_freqs_waits(:,1), hi_tone_freqs), :);
    [r,p] = fit_line(trial_corrs_hi, wait_times_hi);
    for idot = hi_tone_nums
        tone_idx = trialTones_hi(:,1) == unqfrq41(idot);
        plot(trial_corrs_hi(tone_idx),  trialTones_hi(tone_idx, 2), 'o', 'color', tcolors(idot,:));
    end
    title(['r=' num2str(round(r*100)/100) ';p=' num2str(round(p*1000)/1000)])
end
    ylim([0 60]); set(gca,'TickLength',[0, 0]); box off;

    
% poor and rich output

poor_tone_rvals = [trial_corrs_lo;trial_corrs_hi];
rich_tone_rvals = [trial_corrs_med];

poor_tone_waits = [wait_times_lo;wait_times_hi];
rich_tone_waits = [wait_times_med];
    
    
    
    
    
end

function [trl_mtx, trl_idx, frame_times, traces] = load_trl_mtx(subject, current_session)
    beh_folder = ['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_imaging_hpc\' subject];
    if current_session <=6
        sesh_str = 'preprobe';
        probe_num = num2str(current_session);
        session = get_file_paths_targeted(beh_folder, [sesh_str '_0' probe_num], 'LED');
    elseif ismember(current_session, [7 8])
        sesh_str = 'postprobe';
        probe_num = num2str(current_session-6);
        session = get_file_paths_targeted(beh_folder, [sesh_str '_0' probe_num], 'LED');
    else
        sesh_str = 'mevar';
        sesh_num = num2str(current_session-8);
        session = get_file_paths_targeted(beh_folder, [sesh_str '0' sesh_num], '01.mat');
        
    end
    
    load(session{1}, 'trl_mtx', 'trl_idx', 'frame_times', 'traces')
end





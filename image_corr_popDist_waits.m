function image_corr_popDist_waits(subject, target_session, comparison_sessions, activity_mtx, cell_regist_mtx, session_number_idx, trial_number_idx, time_bin_idx, time_series_event_spacing)
% breaks wait times from current probe into thirds (low tones,
% (mevar-rewarded) middle tones, and high tones) and correlates them with
% the average population distance from the target probe to the two input 
% comparison probes


%% load trl mtx
[trl_mtx, trl_idx] = load_trl_mtx(subject, target_session);



%% compute wait times
hm_cell_trials = unique(trl_idx);
hm_cell_trials = intersect(hm_cell_trials, find(trl_mtx(:,3)==0));
probe_freqs_waits = trl_mtx(hm_cell_trials, [2 12]);
probe_freqs_waits(:,1) = floor(probe_freqs_waits(:,1));



%% compute distances
all_mean_diff_distances = nan(size(probe_freqs_waits,1), 1);
for itrl = 1:size(probe_freqs_waits,1)
    [~, all_mean_diff_distances(itrl)] = ...
        image_timeBin_popDist_line(target_session, comparison_sessions, itrl, activity_mtx, cell_regist_mtx, session_number_idx, trial_number_idx, time_bin_idx, time_series_event_spacing); 
    %close % CLOSE FIGURE
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
pop_dists = all_mean_diff_distances(ismember(probe_freqs_waits(:,1), lo_tone_freqs))
wait_times = probe_freqs_waits(ismember(probe_freqs_waits(:,1), lo_tone_freqs), 2)
[r,p] = fit_line(pop_dists, wait_times);
for idot = lo_tone_nums
    plot(all_mean_diff_distances(probe_freqs_waits(:,1) == unqfrq41(idot)),  probe_freqs_waits(probe_freqs_waits(:,1) == unqfrq41(idot), 2), 'o', 'color', tcolors(idot,:));
end
ylim([0 60])
title(['r=' num2str(round(r*100)/100) ';p=' num2str(round(p*1000)/1000)])

% low
subplot(1,3,2); hold on
pop_dists = all_mean_diff_distances(ismember(probe_freqs_waits(:,1), med_tone_freqs));
wait_times = probe_freqs_waits(ismember(probe_freqs_waits(:,1), med_tone_freqs), 2);
[r,p] = fit_line(pop_dists, wait_times);
for idot = med_tone_nums
    plot(all_mean_diff_distances(probe_freqs_waits(:,1) == unqfrq41(idot)),  probe_freqs_waits(probe_freqs_waits(:,1) == unqfrq41(idot), 2), 'o', 'color', tcolors(idot,:));
end
ylim([0 60])
title(['r=' num2str(round(r*100)/100) ';p=' num2str(round(p*1000)/1000)])

% low
subplot(1,3,3); hold on
pop_dists = all_mean_diff_distances(ismember(probe_freqs_waits(:,1), hi_tone_freqs));
wait_times = probe_freqs_waits(ismember(probe_freqs_waits(:,1), hi_tone_freqs), 2);
[r,p] = fit_line(pop_dists, wait_times);
for idot = hi_tone_nums
    plot(all_mean_diff_distances(probe_freqs_waits(:,1) == unqfrq41(idot)),  probe_freqs_waits(probe_freqs_waits(:,1) == unqfrq41(idot), 2), 'o', 'color', tcolors(idot,:));
end
ylim([0 60])
title(['r=' num2str(round(r*100)/100) ';p=' num2str(round(p*1000)/1000)])


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





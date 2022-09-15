function [all_trial_nums, all_freqs, all_waits, all_rich_idx, all_TTL_idx] = ALL_wait_times_sesh_plot(file_keywords)
%plots the mean wait times over session, annotated for trial type

% get all session files
fp = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_optoExp_hpc\';
%fp = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar\';
session_files = get_file_paths_targeted(fp, file_keywords);

selected_problems = [4 5 6];
session_files_hold = [];
for iprob = selected_problems
    session_files_hold = [session_files_hold; session_files(contains(session_files,['mevar0' num2str(iprob)]))];
end
session_files = unique(session_files_hold)


% soft prep
all_trial_nums = [];
all_freqs = [];
all_waits = [];
all_rich_idx = [];
all_TTL_idx = [];

% plot (will happen either way; close later if desired)
figure; hold on

% iterate through session files
for isf = 1:size(session_files,1)
    
    % load session file
    session_files{isf};
    load(session_files{isf}, 'trl_mtx');
    
    % compute
    [trial_nums, freqs, wait_times, rich_idx, TTL_idx] = wait_times_sesh_plot(trl_mtx);
    
    % normalize by prior probe wait times
    prob_num = session_files{isf}(strfind(session_files{isf}, 'mevar0')+6);
    
    % load
    all_trial_nums = [all_trial_nums; trial_nums];
    all_freqs = [all_freqs; freqs];
    all_waits = [all_waits; wait_times];
    all_rich_idx = [all_rich_idx; rich_idx];
    all_TTL_idx = [all_TTL_idx; TTL_idx];

end

% line plots
colors = [.4 .4 .4; 0.05 0.29 0.65];
figure; hold on

    %legend
    legend_trick(colors, '-')
    legend_trick(colors, '--')

idx = all_rich_idx==1 & all_TTL_idx==0; dots_to_line_plot(all_trial_nums(idx), all_waits(idx), 21, colors(1,:), '-') % RichTone LaserOFF
idx = all_rich_idx==1 & all_TTL_idx==1; dots_to_line_plot(all_trial_nums(idx), all_waits(idx), 21, colors(2,:), '-') % RichTone LaserON
idx = all_rich_idx==0 & all_TTL_idx==0; dots_to_line_plot(all_trial_nums(idx), all_waits(idx), 21, colors(1,:), '--') % PoorTone LaserOFF
idx = all_rich_idx==0 & all_TTL_idx==1; dots_to_line_plot(all_trial_nums(idx), all_waits(idx), 21, colors(2,:), '--') % PoorTone LaserON

% aesthetics
set(gca,'TickLength',[0, 0]);
ylabel('Wait Durations (s)')
xlabel('Trial number')
xlim([-2.01 125.01])
legend({'RichTone LaserOFF','RichTone LaserON','PoorTone LaserOFF','PoorTone LaserON'}, 'location', 'northeastoutside')
ylim([0 60])

end

function dots_to_line_plot(trial_nums, waits, binsize, color, line_type)

% preallocate
bin_means = [];
bin_ses = [];
xaxis_tn = [];

% compute
count = 0;
for itrl = 1 : 1 : max(trial_nums)
    count = count+1;
    
    binrng = binsize/2;
    
    local_waits = waits(trial_nums>itrl-binrng & trial_nums<itrl+binrng);
    bin_means = [bin_means; mean(local_waits)];
    bin_ses = [bin_ses; std(local_waits)./sqrt(sum(~isnan(local_waits)))];
    
end

% smooth
%bs = 3;
%bin_means(~isnan(bin_means)) = smooth(bin_means(~isnan(bin_means)), bs);
%bin_ses(~isnan(bin_ses)) = smooth(bin_ses(~isnan(bin_ses)), bs);

% plot
mtm = max(trial_nums);
plot(1 : mtm, bin_means, line_type, 'linewidth', 2, 'color', color)
plot(1 : mtm, bin_means+bin_ses, line_type, 'linewidth', 1, 'color', color)
plot(1 : mtm, bin_means-bin_ses, line_type, 'linewidth', 1, 'color', color)

xticks([0 mtm/4 (mtm/4)*2 (mtm/4)*3 mtm]); xticklabels(xticks/mtm); xlabel('Proportion of session')

end
function all_waits = plot_discrim_curves(datafolder, subject, varargin)
%input string matching the subject file and additional string info to
%identify the desired trials (desired files should contain string).

% folderpath
folderpath = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\';
folderpath = [folderpath datafolder];
folderpath = [folderpath '\' subject];

% filter with additional input string constraints
fpaths = get_file_paths_targeted(folderpath, varargin)
if isempty(fpaths)
    error('No matching files found')
end

% colors
%colors = distinguishable_colors(length(fpaths));
colors = winter(length(fpaths));

% preallocate waits
all_waits = cell(length(fpaths),3);
mean_waits = nan(length(fpaths),3);
se_waits = nan(length(fpaths),3);
zdiffs = nan(length(fpaths),1);

% legend prep
figure; hold on
legend_input = cell(length(fpaths),1);
for ifile = 1:length(fpaths)
   plot(6500,1,'-', 'color', colors(ifile,:), 'linewidth', 1.5)
   legend_input{ifile} = ['d' fpaths{ifile}(end-5:end-4)];
end

% iteratively load files
unq_frq_sesh = [];
for ifile = 1:length(fpaths)
    load(fpaths{ifile}, 'trl_mtx', 'medass_cell');
    
    %zscore waits
    %trl_mtx(trl_mtx(:,3)==0,12) = zscore_mtx(trl_mtx(trl_mtx(:,3)==0,12));
    
    % compute waits
    [mwd, mwd_freqs] = wait_times_prep(trl_mtx, 1); %1 for means, 2 for all
    [wd, wd_freqs] = wait_times_prep(trl_mtx, 2); %1 for means, 2 for all
    unq_frq = unique(wd_freqs);
    unq_frq_sesh = [unq_frq_sesh; unq_frq'];
    
    % reward probability distribution
    [prob_dist, pd_freq] = rwd_prob_by_freq(medass_cell);
    rich_tones = pd_freq(prob_dist>.5);
    
    %z differences
    
    if ~isempty(wd(ismember(wd_freqs,rich_tones))) && ~isempty(wd(~ismember(wd_freqs,rich_tones)))
        zdiffs(ifile) = zdiff(wd(ismember(wd_freqs,rich_tones)), wd(~ismember(wd_freqs,rich_tones)));
    end
    % resort wait output
    %[unq_frq,frq_sort_idx] = sort(unq_frq);
    %mwd = mwd(frq_sort_idx);
    %wd = wd(frq_sort_idx);
    
    % load
    mean_waits(ifile,1:length(mwd)) = mwd;
    for ifrq = 1:length(unq_frq)
        current_frq = unq_frq(ifrq);
        all_waits{ifile,ifrq} = wd(wd_freqs==current_frq);
        se_waits(ifile,ifrq) = std(wd(wd_freqs==current_frq))./sqrt(sum(~isnan(wd(wd_freqs==current_frq))));
    end

end

% overall unqfrq
unq_frq = unique(unq_frq_sesh);

% plot all waits
for ifile = 1:size(mean_waits,1)
    
    % individual waits
    for ifrq = 1:size(unq_frq_sesh,2)

        % jitter x positions
        jf_hold = jitter_xpos(repmat(unq_frq_sesh(ifile,ifrq), size(all_waits{ifile,ifrq})), all_waits{ifile,ifrq}, unq_frq_sesh(ifile,ifrq).*0.50);
        
        % plot
        plot(jf_hold, all_waits{ifile,ifrq}, 'o', 'color', colors(ifile,:), 'markersize', 7);
        
    end
end

% plot all curves
for ifile = 1:size(mean_waits,1)
        errorbar(unq_frq_sesh(ifile,:)+((rand(1)-.5)*0.125*unq_frq_sesh(ifile,:)), mean_waits(ifile,1:length(mwd)), se_waits(ifile,1:length(mwd)), '-', 'color', colors(ifile,:), 'linewidth', 1.5)
end

% aesthetics
set(gca,'TickLength',[0, 0]);
xlabel('Tone Frequency (Hz)')
ylabel('Wait Durations (s)')
set(gca, 'XScale', 'log')
title('Discrimination curves')
ylim_hold = ylim;
if ylim_hold(1)>0
    ylim_hold(1) = 0;
end
if ylim_hold(2)>60
     ylim_hold(2) = 60;
end
ylim(ylim_hold)
xlim([4500 42000])
xticks([5000 8500 14000 23000 35000])
legend(legend_input, 'location', 'northeastoutside')

% zdiff learning curve
figure

zdiffs

plot(zdiffs, '-o')
xlim([0.5 length(zdiffs)+0.5])
set(gca,'TickLength',[0, 0]); box off;
hold on; plot(xlim, [1 1].*0, 'k--')
xticks(1:length(zdiffs)); xlabel('Day')
ylabel('Preference for rich tone (zdiff)')

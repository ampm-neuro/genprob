function plot_all_waits(fpaths, varargin)
% paths to trl_mtx files
% can also input a plot color

% check input
if nargin==2
    colors = varargin{1};
end

% iterate through fpaths
all_wait_times = [];
all_frequencies = [];
for ipath = 1:size(fpaths,1)
        
    % load
    load(fpaths{ipath}, 'trl_mtx')

    % compute waits
    means_or_all = 2; %means=1, all=2
    [wait_times, frequencies] = wait_times_prep(trl_mtx, 1);
    
    % plot
    hold on
    if exist('colors','var')
        plot(frequencies, wait_times, 'o', 'color', colors)
    else
        plot(frequencies, wait_times, 'o', 'color', [1 1 1].*0.5)
    end
    
    % load
    all_wait_times = [all_wait_times; wait_times];
    all_frequencies = [all_frequencies; frequencies];

end


% plot errorbar
unq_freq = unique(all_frequencies);
all_means = nan(length(unq_freq),1);
all_ses = nan(length(unq_freq),1);
for itone = 1:length(unq_freq)
   relv_freq = unq_freq(itone);
   freq_idx = all_frequencies == relv_freq;
   tone_waits = all_wait_times(freq_idx);
   all_means(itone) = nanmean(tone_waits);
   all_ses(itone) = nanstd(tone_waits)./sum(~isnan(tone_waits));
end

if exist('colors','var')
    errorbar(unq_freq, all_means, all_ses, '-', 'color', colors, 'linewidth', 2)
    %plot(unq_freq, smooth(all_means), '-', 'color', colors, 'linewidth', 2)
else
    plot(unq_freq, smooth(all_means), '-', 'color', [1 1 1].*0.5, 'linewidth', 2)
end



% aesthetics
set(gca,'TickLength',[0, 0]);
xlabel('Tone Frequency (Hz)')
ylabel('Wait Durations (s)')
set(gca, 'XScale', 'log')
xlim([4500 42000])
xticks([5000 8500 14000 23000 35000])
ylim([0 60])
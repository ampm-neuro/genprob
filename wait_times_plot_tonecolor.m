function [h1, unq_frq, mean_wait_times] = wait_times_plot_tonecolor(trl_mtx)
%plots the mean wait times at each frequency
plot_type=1;


% rich tones
load('unqfrq41.mat', 'unqfrq41');
tcolors = parula(71); tcolors = tcolors(26:66,:);

% all wait times and corresponding frequencies
if length(unique(floor(trl_mtx(~isnan(trl_mtx(:,2)),2))))>10
    [wait_times, frequencies] = wait_times_prep(trl_mtx, 2, 0); %probe
else
    [wait_times, frequencies] = wait_times_prep(trl_mtx, 2, 0);
end

% means and ses
unq_frq = unique(frequencies);
mean_wait_times = nan(length(unq_frq),1);
se_wait_times = nan(size(mean_wait_times));
wt = cell(1,length(unq_frq));
for ifrq = 1:length(unq_frq)
    wt{ifrq} = wait_times(frequencies==unq_frq(ifrq));
    mean_wait_times(ifrq) = mean(wait_times(frequencies==unq_frq(ifrq)));
    se_wait_times(ifrq) = std(wait_times(frequencies==unq_frq(ifrq)))/sqrt(sum(frequencies==unq_frq(ifrq)));
end

% stats
%{
if length(unq_frq) == 2
    [a b c d] = ttest2(wt{1}, wt{2})
end
%}


% plot individual waits
if ismember(plot_type, [1 3])
    hold on; 
    for ifrq = 1:length(unq_frq)
            if length(wait_times(frequencies==unq_frq(ifrq)))>1
                xpos_jit = jitter_xpos(repmat(unq_frq(ifrq), size(wait_times(frequencies==unq_frq(ifrq)))),wait_times(frequencies==unq_frq(ifrq)), unq_frq(ifrq)/2);
            else
                xpos_jit = repmat(unq_frq(ifrq), size(wait_times(frequencies==unq_frq(ifrq))));
            end
        if ifrq==1
            h1 = plot(xpos_jit, wait_times(frequencies==unq_frq(ifrq)), 'o', 'color', tcolors(ismember(unqfrq41, unq_frq(ifrq)),:));
        else
            plot(xpos_jit, wait_times(frequencies==unq_frq(ifrq)), 'o', 'color', tcolors(ismember(unqfrq41, unq_frq(ifrq)),:));
        end
    end
end

% plot means and se
hold on

if ismember(plot_type, [2 3])
    errorbar(unq_frq, mean_wait_times, se_wait_times, 'linewidth', 1.5, 'color', tcolors(ismember(unqfrq41, unq_frq(ifrq)),:));
end
    

% aesthetics
set(gca,'TickLength',[0, 0]);
xlabel('Tone Frequency (Hz)')
ylabel('Wait Durations (s)')
ylim_hold = ylim; ylim([0 ylim_hold(2)]);
set(gca, 'XScale', 'log')
title('Waits')
xlim([4500 42000])
xticks([5000 8500 14000 23000 35000])

    




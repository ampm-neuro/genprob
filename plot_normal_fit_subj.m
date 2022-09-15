function [coefEsts, fit_yvals, obs_yvals] = plot_normal_fit_subj(trl_mtx, varargin)
% overlays normal distribution on single subject wait times

if nargin == 2 
    plot_on = varargin{1};
elseif nargin == 3
    plot_on = varargin{1};
    colors = varargin{2};
else
    plot_on = 3;
end

%should probably use means because otherwise the freqs with more wait
%samples will disproportionately bias the fit
means_or_all = 1; % 1 = means, 2 = all waits

% get wait times
[wait_times, frqs, frq_nums] = wait_times_prep(trl_mtx, means_or_all, 0);
obs_yvals = wait_times;

%zscore waits
%waits = zscore_mtx(waits);

% fit waits
if ismember(plot_on, [1 3]) && exist('colors', 'var')
    hold on
    h1 = plot(frqs, wait_times, 'o', 'color', colors);
elseif ismember(plot_on, [1 3])
    hold on
    h1 = plot(frqs, wait_times, 'o');
end

% fit normal distribution
[fit_yvals, coefEsts, modelFun] = ampm_normal_fit(frq_nums, wait_times);

% interp a smooth curve
unqfrq = csvread('unq_frq_41.csv');
unqfrq_interp = interp1(1:length(unqfrq), unqfrq, (1:0.001:length(unqfrq))', 'spline');
fitys_interp = modelFun(coefEsts, (1:0.001:length(unqfrq))');

% plot 
if ismember(plot_on, [1 2 3])

    % overlay 
    if plot_on==3 && exist('colors', 'var')
        %line(unique(frqs), fit_yvals, 'color', h1.Color)
        line(unqfrq_interp, fitys_interp, 'color', h1.Color)
    elseif plot_on==2 && exist('colors', 'var')
        %line(unique(frqs), fit_yvals, 'color', colors)
        line(unqfrq_interp, fitys_interp, 'color', colors)
    elseif ismember(plot_on, [2 3])
        %line(unique(frqs), fit_yvals)
        line(unqfrq_interp, fitys_interp)
    end
    set(gca, 'XScale', 'log')
    set(gca,'TickLength',[0, 0]); box off;

    %aesthetics
    set(gca,'TickLength',[0, 0]); 
    set(gca, 'XScale', 'log');
    box off;
    ylim_hold = ylim;
    ylim([0 ylim_hold(2)])
    xlim([4650 38000])
    xticks([5000 8500 14000 23000 35000])
    xticklabels({'5000', '8500', '14000', '23000', '35000'})
    xlabel('Tone Frequency (Hz)')
    ylabel('Mean wait times (s)')
    %rich_bounds; %script plotting lines around rich area

end




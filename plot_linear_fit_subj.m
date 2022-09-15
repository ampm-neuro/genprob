function [intercept, slope, fit_yvals, obs_yvals, gof] = plot_linear_fit_subj(wd, varargin)
% overlays normal distribution on single subject wait times

if nargin == 2 
    plot_on = varargin{1};
elseif nargin == 3
    plot_on = varargin{1};
    colors = varargin{2};
else
    plot_on = 1;
end

%should probably use means because otherwise the freqs with more wait
%samples will disproportionately bias the fit
means_or_all = 1; % 1 = means, 2 = all waits

%hard code frequencies 
unq_frq = [5000,5249,5511,5786,6074,6377,6695,7028,7379,7747,8133,8538,...
    8964,9411,9880,10372,10890,11432,12002,12601,13229,13888,14581,...
    15307,16070,16872,17713,18596,19523,20496,21518,22590,23716,24899,...
    26140,27443,28811,30247,31755,33338,35000];

% convert inputs into vectors
if means_or_all == 1 %mean waits only
    waits = cellfun(@mean, wd);
    frqs = []; 
    frq_nums = [];
    for ifrq = 1:length(unq_frq)
        frqs = [frqs; repmat(unq_frq(ifrq),size(mean(wd{ifrq})))];
        frq_nums = [frq_nums; repmat(ifrq,size(mean(wd{ifrq})))];
    end

elseif means_or_all == 2 %all waits
    waits = cell2mat(wd);
    frqs = []; 
    frq_nums = [];
    for ifrq = 1:length(unq_frq)
        frqs = [frqs; repmat(unq_frq(ifrq),size(wd{ifrq}))];
        frq_nums = [frq_nums; repmat(ifrq,size(wd{ifrq}))];
    end
else
    error('must hard code to use mean waits or all waits')
end


%zscore waits
%waits = zscore_mtx(waits);

% fit waits

if ismember(plot_on, [1 3]) && exist('colors', 'var')
    hold on
    h1 = plot(frqs, waits, 'o', 'color', colors);
elseif ismember(plot_on, [1 3])
    hold on
    h1 = plot(frqs, waits, 'o');
end

% fit linear distribution
obs_yvals = waits;
model = fitlm(frq_nums,obs_yvals,'linear');
intercept = model.Coefficients.Estimate(1);
slope = model.Coefficients.Estimate(2);
slope_pval = model.Coefficients.pValue(2);

fit_yvals = intercept + frq_nums.*slope;

% goodness of fit
gof = goodnessOF(obs_yvals, fit_yvals);

% plot 
if ismember(plot_on, [1 2 3])

    % overlay 
    if plot_on==3 && exist('colors', 'var')
        line(unq_frq, fit_yvals, 'color', h1.Color)
    elseif plot_on==2 && exist('colors', 'var')
        line(unq_frq, fit_yvals, 'color', colors)
    elseif ismember(plot_on, [2 3])
        line(unq_frq, fit_yvals)
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
    rich_bounds; %script plotting red lines around rich area

end




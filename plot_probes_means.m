function [fixed_effects, random_effects, beta_pvals, all_subjects] = plot_probes_means(datafolder, probe_num, day_delay, new_learn_num, varargin)
% plots curve using subject means at each tone frequency
% probe_num is 1 2 and/or 3 (first, second, or third probe session)
% delay day is 1 and/or 30 (delay between last training day and first
% probe)

%input checks
if size(probe_num,1)>size(probe_num,2)
    probe_num = probe_num';
end
if size(day_delay,1)>size(day_delay,2)
    day_delay = day_delay';
end
if size(new_learn_num,1)>size(new_learn_num,2)
    new_learn_num = new_learn_num';
end

if nargin ==5
    plot_on = varargin{1};
elseif nargin == 6
    plot_on = varargin{1};
    colors = varargin{2};
else
    plot_on = 1; %default is plot ON
end

% folderpath
folderpath = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\';
folderpath = [folderpath datafolder]


%get probe files
new_learning_strings = [];
for i = new_learn_num
    new_learning_strings = [new_learning_strings {['learn0' num2str(i)]}];
end

[probe_files_constrained, all_subjects] = get_file_paths_targeted_II(folderpath, {'probe', ['_0' num2str(probe_num) '_']}, {[num2str(day_delay) 'd'], new_learning_strings});
if isempty(probe_files_constrained)
    return
end

% plot probe files
[fixed_effects, random_effects, beta_pvals, ~, modelFun, unqfrq] = mixed_model_fit(probe_files_constrained, 0);


%ses of all random_effects
ses_RE = nanstd(random_effects,[],2);
ses_RE = ses_RE./sqrt(sum(~isnan(random_effects)));


if plot_on > 0
    hold on
    
    %plot mean and se of all probe fit curves
    if exist('colors', 'var')
        h1 = plot(unqfrq, modelFun(fixed_effects, 1:length(unqfrq)), '-', 'linewidth', 2, 'color', colors);
    else
        h1 = plot(unqfrq, modelFun(fixed_effects, 1:length(unqfrq)), '-', 'linewidth', 2);
    end
    try
        plot(unqfrq, modelFun(fixed_effects + ses_RE, 1:length(unqfrq)), '-', 'linewidth', 1, 'color', h1.Color)
        plot(unqfrq, modelFun(fixed_effects - ses_RE, 1:length(unqfrq)), '-', 'linewidth', 1, 'color', h1.Color)
    catch
    end
    
    set(gca,'TickLength',[0, 0]); 
    set(gca, 'XScale', 'log');
    box off;
    ylim auto
    %ylim_hold = ylim;
    %ylim([0 ylim_hold(2)])
    xlim([4650 38000])
    xticks([5000 8500 14000 23000 35000])
    xticklabels({'5000', '8500', '14000', '23000', '35000'})
    xlabel('Tone Frequency (Hz)')
    ylabel('Mean wait times (s)')
    rich_bounds; %script plotting red lines around rich area
    legend({'mean fit', 'std error'}, 'location', 'northeast')

end






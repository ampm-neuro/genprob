function [fixed_effects, random_effects, beta_pvals, coef_names, modelFun, unqfrq] = mixed_model_fit_normal_logistic(probe_paths, varargin)
% uses probe data from all subjects and to fit a nonlinear mitonesed-effects...
% regression model and returns estimates of the fitonesed effects
%
% plot on :
% 1 for mean line, 
% 2 for mean line + subj mean  lines, 
% 3 for mean line + subj mean lines + wait times scatter

% check number of paths
if isempty(probe_paths)
    disp('no paths')
    return
end

% plot on
if nargin ==2
    plot_on = varargin{1}; 
elseif nargin ==3
    plot_on = varargin{1}; 
    colors = varargin{2};
else
    plot_on = 0;
end

% number of subjects
num_subjs = size(probe_paths,1);
session_ct_idx = 1:num_subjs;

% colors
if ~exist('colors', 'var')
    colors = distinguishable_colors(num_subjs);
end

% compute all wait times
all_wait_times = cell(num_subjs,1);
for isesh = session_ct_idx
    load(probe_paths{isesh}, 'trl_mtx')
    [trial_wait_times, trial_freqs] = wait_times_prep(trl_mtx, 1, 2); %1 for means, 2 for all
    all_wait_times{isesh} = [trial_wait_times, trial_freqs];
end

% unique tones
unqfrq = csvread('unq_frq_41.csv');
all_tones = (1:length(unqfrq))';

% specify model
%{
modelFun = @(coefEsts,tones) coefEsts(4).*exp(-((tones-coefEsts(5))./coefEsts(6)).^2)...
    + (coefEsts(2).*(exp(tones-coefEsts(3))./(exp(tones-coefEsts(3))+1)))...
    + coefEsts(1);
%}
modelFun = @(coefEsts,tones) coefEsts(4).*exp(-((tones-coefEsts(3))./coefEsts(4)).^2)...
    + (coefEsts(2).*(exp(tones-coefEsts(3))./(exp(tones-coefEsts(3))+1)))...
    + coefEsts(1);

% prepare inputs
wait_input = [];
subject_input = [];
subj_cell = cell(num_subjs,1);
tone_input = [];
for isesh = session_ct_idx
    
    unq_tones = unique(all_wait_times{isesh}(:,2));
    numsamps = histc(all_wait_times{isesh}(:,2), unq_tones);
    for itone = 1:length(unique(all_wait_times{isesh}(:,2)))
        subj_cell{isesh} = [subj_cell{isesh}; repmat(isesh, numsamps(itone), 1)];
    end    
    
    subject_input = [subject_input; subj_cell{isesh}];
    tone_input = [tone_input; all_wait_times{isesh}(:,2)];
    wait_input = [wait_input; all_wait_times{isesh}(:,1)];
        
end

% set tones to numbers
[~,~,comb_tones] = unique([tone_input; unqfrq']);
tone_input = comb_tones(1:length(tone_input));

% remove nans
nnan_idtones = ~isnan(tone_input) & ~isnan(wait_input) & ~isnan(subject_input);
tone_input = tone_input(nnan_idtones);
wait_input = wait_input(nnan_idtones);
subject_input = subject_input(nnan_idtones);

% mitonesed effects model
%starting_vals = [1 1 1 1 1 1];
%starting_vals = [12.6 3 20.5 10 20.5 9];
starting_vals = [12.6 3 20.5 10];


[fixed_effects,~,stats,random_effects] = nlmefit(tone_input, wait_input, subject_input, [], modelFun, starting_vals, 'RefineBeta0', 'on');

% output     
coef_names = {'intercept', 'amplitude', 'mu', 'sigma'};
beta_pvals = 1 - tcdf( abs(fixed_effects')./stats.sebeta , stats.dfe);

% plot model
%

if plot_on > 0
    figure; hold on
end

    % plot each subject's wait times and model
    isubj_ct = 0;
    for isubj = session_ct_idx
        isubj_ct = isubj_ct + 1;

        % subject wait times
        if plot_on == 3
            plot(unqfrq(tone_input(subject_input==isubj_ct)), all_wait_times{isubj}, 'o', 'color', colors(isubj,:))
        end

        % subject model
        if ismember(plot_on, [2 3])
            PHI = fixed_effects + random_effects(:,isubj_ct);
            %; %fixed effects + random effects
            plot(unqfrq,modelFun(PHI,all_tones), 'color', colors(isubj,:))
        end

    end
    
    % plot full model
    if ismember(plot_on, [1 2 3])
        if nargin==3
            plot(unqfrq, modelFun(fixed_effects, all_tones), '-', 'color', colors(end,:), 'linewidth', 3)
        else
            plot(unqfrq, modelFun(fixed_effects, all_tones), 'k-', 'linewidth', 3)
        end
    
        % aesthetics
        set(gca, 'XScale', 'log')
        xlim([4600 38000])
        ylim([0 40])
        set(gca,'TickLength',[0, 0]); box off;
        %rich_bounds
    end
    
    %{
    figure; hold on
    for i = 1:size(colors,1)
        plot(i,1,'.','color', colors(i,:), 'markersize', 30)
    end
    %}




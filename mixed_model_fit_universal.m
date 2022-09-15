function [fixed_effects, random_effects, beta_pvals, coef_names, modelFun, unqfrq] = mixed_model_fit_universal(probe_paths, varargin)
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
else
    plot_on = 0;
end

% number of sessions
num_paths = size(probe_paths,1);

% included sessions
session_ct_idx = 1:num_paths;
probe_paths = probe_paths(session_ct_idx);

% Array of covariate values
    d_idx = strfind(probe_paths{1},'d');
    d_idx = d_idx(end);

    delays = [];
    for i = 1:length(probe_paths)
        delays = [delays; str2double(probe_paths{i}(d_idx-5 : d_idx-4))];
    end
    delays = delays -1 ;
    pos_delays = 0:max(delays);

% further constrain sessions
    %delays = delays(session_ct_idx);
    %delays(delays==1) = 0;
    %delays(delays==30) = 1;
    %session_ct_idx = ismember(delays, [0:6]);
    %probe_paths = probe_paths(session_ct_idx);
    %num_paths = size(probe_paths,1);
    %delays = delays(session_ct_idx);

    %probe_paths
    probe_paths

    
% get unique subject ids
unq_subj = [];
for ipath = 1:size(probe_paths,1)
    gen_loc = strfind(probe_paths{ipath},'gen'); gen_loc = max(gen_loc);
    unq_subj = [unq_subj; probe_paths{ipath}(gen_loc-9 : gen_loc-2)];
end
[~, ~, unq_subj] = unique(unq_subj, 'rows', 'stable');

% colors
colors = distinguishable_colors(num_paths);

% compute all wait times
all_wait_times = cell(num_paths,1);
for isesh = 1:num_paths
    load(probe_paths{isesh}, 'trl_mtx','medass_cell')
    [trial_wait_times, trial_freqs] = wait_times_prep(trl_mtx, 2); %1 for means, 2 for all
    all_wait_times{isesh} = [trial_wait_times, trial_freqs];
end

% unique tones
unqfrq = csvread('unq_frq_41.csv');
all_tones = (1:0.001:length(unqfrq))';
unqfrq_interp = interp1(1:length(unqfrq), unqfrq, all_tones, 'spline');

% specify model
%modelFun = @(coefEsts,tones) coefEsts(2).*exp(-((tones-coefEsts(3))./coefEsts(4)).^2) + coefEsts(1);
modelFun = @(coefEsts,tones,delay,subject) (coefEsts(2) + delay*coefEsts(6)).*exp(-((tones-(coefEsts(3) + delay*coefEsts(7)))./(coefEsts(4) + delay*coefEsts(8))).^2) + (coefEsts(1) + delay*coefEsts(5));

% prepare inputs
wait_input = [];
sesh_input = [];
sesh_cell = cell(num_paths,1);
tone_input = [];
tone_cell = cell(num_paths,1);
delay_input = [];
delay_cell = cell(num_paths,1);
for isesh = 1:num_paths
    
    % delay condition
    after_delay_cond = strfind(probe_paths{isesh},'d-');
    delay_condition = str2double(probe_paths{isesh}(after_delay_cond-2 : after_delay_cond-1));
    
    % tones and subjects
    unq_tones = unique(all_wait_times{isesh}(:,2));
    numsamps = histc(all_wait_times{isesh}(:,2), unq_tones);
    for itone = 1:length(unique(all_wait_times{isesh}(:,2)))
        tone_cell{isesh} = [tone_cell{isesh}; repmat(itone, numsamps(itone), 1)];
        sesh_cell{isesh} = [sesh_cell{isesh}; repmat(isesh, numsamps(itone), 1)];
        delay_cell{isesh} = [delay_cell{isesh}; repmat(delay_condition, numsamps(itone), 1)];
    end    
    
    sesh_input = [sesh_input; sesh_cell{isesh}];
    tone_input = [tone_input; all_wait_times{isesh}(:,2)];
    wait_input = [wait_input; all_wait_times{isesh}(:,1)];
    delay_input = [delay_input; delay_cell{isesh}];
        
end

% set tones to numbers
[~,~,comb_tones] = unique([tone_input; unqfrq']);
tone_input = comb_tones(1:length(tone_input));
    
% remove nans
nnan_idtones = ~isnan(tone_input) & ~isnan(wait_input) & ~isnan(sesh_input);
tone_input = tone_input(nnan_idtones);
wait_input = wait_input(nnan_idtones);
sesh_input = sesh_input(nnan_idtones);
delay_input = delay_input(nnan_idtones);

% starting values
%starting_vals = [1 1 1 1 0 0 0 0];
starting_vals = [12.6 10 15 9 0 0 0 0];

% model
[fixed_effects,re_var,stats,random_effects] = nlmefit(tone_input,wait_input,sesh_input,delays,modelFun,starting_vals);
re_var

% output     
coef_names = {'intercept', 'amplitude', 'mu', 'sigma'};
%beta_pvals = 0;
beta_pvals = 1 - tcdf( abs(fixed_effects')./stats.sebeta , stats.dfe);


% plot model
%

% colors
if ~exist( 'pos_delays', 'var')
    pos_delays = [0 1]
end
colors = distinguishable_colors(length(pos_delays));


% figure on
figure; hold on
legend_trick(colors, '-')

    % plot each subject's wait times and model
    isubj_ct = 0;
    for sesh = 1:length(delays)
        isubj_ct = isubj_ct + 1;

        % subject wait times
        if plot_on == 3
              % plot(unqfrq(tone_cell{isubj}), cell2mat(all_wait_times{isubj}), 'o', 'color', colors(isubj,:))
              plot(unqfrq(tone_cell{sesh}), all_wait_times{sesh}(:,1), 'o', 'color', colors(delays(isubj_ct)+1,:))

        end

        % subject model
        if ismember(plot_on, [2 3])
            PHI = fixed_effects + random_effects(:,isubj_ct);
            %plot(unqfrq,modelFun(PHI,all_tones, delays(isubj_ct)), 'color', colors(isubj,:))
                plot(unqfrq_interp,modelFun(PHI, all_tones, delays(isubj_ct)), 'color', colors(delays(isubj_ct)+1,:))
        end

    end
    
    % plot full model
    if ismember(plot_on, [1 2 3])       
        plot(unqfrq_interp, modelFun(fixed_effects, all_tones, unique(delays)'), 'k-', 'linewidth', 3)
    
        % aesthetics
        set(gca, 'XScale', 'log')
        xlim([4600 38000])
        ylim([0 60])
        set(gca,'TickLength',[0, 0]); box off;
        rich_bounds
    end
    
    figure; hold on
    legend_trick(colors, '-')

    if ismember(plot_on, [1 2 3])    
        unq_delays = unique(delays);
        del_ct = 0;
        for unq_del = unique(delays)'
            del_ct = del_ct+1;
            
            plot(unqfrq_interp, modelFun(fixed_effects, all_tones, unq_del), '-', 'color', colors(pos_delays==unq_delays(del_ct),:), 'linewidth', 4)
            plot(unqfrq_interp, modelFun(fixed_effects + std(random_effects(:,delays==unq_del),[],2)./sqrt(sum(delays==unq_del)), all_tones, unq_del), '-', 'color', colors(pos_delays==unq_delays(del_ct),:), 'linewidth', 2)
            plot(unqfrq_interp, modelFun(fixed_effects - std(random_effects(:,delays==unq_del),[],2)./sqrt(sum(delays==unq_del)), all_tones, unq_del), '-', 'color', colors(pos_delays==unq_delays(del_ct),:), 'linewidth', 2)
        end
        % aesthetics
        set(gca, 'XScale', 'log')
        xlim([4600 38000])
        ylim([0 60])
        set(gca,'TickLength',[0, 0]); box off;
        rich_bounds
        legend({'1', '2', '3', '4', '5', '6', '7'})
    end




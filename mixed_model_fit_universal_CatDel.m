function [fixed_effects, random_effects, beta_pvals, coef_names, modelFun, unqfrq] = mixed_model_fit_universal_CatDel(probe_paths, varargin)
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
    session_ct_idx = [1:num_paths];
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
%{
    delays = delays(session_ct_idx);
    session_ct_idx = ismember(delays, [0:6]);
    probe_paths = probe_paths(session_ct_idx);
    num_paths = size(probe_paths,1);
    delays = delays(session_ct_idx);
%}


% print paths
probe_paths

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
%modelFun = @(coefEsts,tones,delay,subject) (coefEsts(2) + delay*coefEsts(6)).*exp(-((tones-(coefEsts(3) + delay*coefEsts(7)))./(coefEsts(4) + delay*coefEsts(8))).^2) + (coefEsts(1) + delay*coefEsts(5));

% delay as a categorical var
%
modelFun = @(coefEsts,tones,group_mtx) sum([coefEsts(2), group_mtx.*coefEsts(6+(size(group_mtx,2)-1)*1:6+(size(group_mtx,2)-1)*2)])...
    .*exp(-((tones-sum([coefEsts(3), group_mtx.*coefEsts(7+(size(group_mtx,2)-1)*2:7+(size(group_mtx,2)-1)*3)]))...
    ./sum([coefEsts(4), group_mtx.*coefEsts(8+(size(group_mtx,2)-1)*3:8+(size(group_mtx,2)-1)*4)])).^2) + sum([coefEsts(1), group_mtx.*coefEsts(5:5+(size(group_mtx,2)-1)*1)]);
    % coefs 1:4 are {'intercept', 'amplitude', 'mu', 'sigma'}. Subsequent
    % coefs are four groups of size(group_mtx,2) coefs describing the
    % relationship between the predictors' 'intercept's (first group), ...
    % 'amplitude's (second group), 'mu's (third group), ...
    % and 'sigma's (last group)} with group 1. 
%}
   
%{
% delay as continuous, subject as categorical
modelFun = @(coefEsts,tones,group_mtx) sum([coefEsts(2), group_mtx(:,1).*coefEsts(6), group_mtx.*coefEsts(10+(size(group_mtx,2)-1)*1:10+(size(group_mtx,2)-1)*2)])...
    .*exp(-((tones-sum([coefEsts(3), group_mtx(:,1).*coefEsts(7), group_mtx.*coefEsts(11+(size(group_mtx,2)-1)*2:11+(size(group_mtx,2)-1)*3)]))...
    ./sum([coefEsts(4), group_mtx(:,1).*coefEsts(8), group_mtx.*coefEsts(12+(size(group_mtx,2)-1)*3:12+(size(group_mtx,2)-1)*4)])).^2) + sum([coefEsts(1), group_mtx(:,1).*coefEsts(5), group_mtx.*coefEsts(9:9+(size(group_mtx,2)-1)*1)]);
    % coefs 1:4 are {'intercept', 'amplitude', 'mu', 'sigma'}. Subsequent
    % coefs  are size(group_mtx,2) groups of 4 coefs describing the
    % relationship between the predictors {'intercept', 'amplitude',...
    % 'mu', 'sigma'} in groups 2:size(gorup_mtx,2) with group 1. 

%}




% prepare inputs
wait_input = [];
sesh_input = [];
sesh_cell = cell(num_paths,1);
tone_input = [];
tone_cell = cell(num_paths,1);
%delay_input = [];
%delay_cell = cell(num_paths,1);
for isesh = 1:num_paths
    
    % delay condition
    after_delay_cond = strfind(probe_paths{isesh},'d-');
    %delay_condition = str2double(probe_paths{isesh}(after_delay_cond-2 : after_delay_cond-1));
    
    % tones and subjects
    unq_tones = unique(all_wait_times{isesh}(:,2));
    numsamps = histc(all_wait_times{isesh}(:,2), unq_tones);
    for itone = 1:length(unique(all_wait_times{isesh}(:,2)))
        tone_cell{isesh} = [tone_cell{isesh}; repmat(itone, numsamps(itone), 1)];
        sesh_cell{isesh} = [sesh_cell{isesh}; repmat(isesh, numsamps(itone), 1)];
        %delay_cell{isesh} = [delay_cell{isesh}; repmat(delay_condition, numsamps(itone), 1)];
    end    
    
    sesh_input = [sesh_input; sesh_cell{isesh}];
    tone_input = [tone_input; all_wait_times{isesh}(:,2)];
    wait_input = [wait_input; all_wait_times{isesh}(:,1)];
    %delay_input = [delay_input; delay_cell{isesh}];
        
end

% set tones to numbers
[~,~,comb_tones] = unique([tone_input; unqfrq']);
tone_input = comb_tones(1:length(tone_input));
    
% remove nans
nnan_idtones = ~isnan(tone_input) & ~isnan(wait_input) & ~isnan(sesh_input);
tone_input = tone_input(nnan_idtones);
wait_input = wait_input(nnan_idtones);
sesh_input = sesh_input(nnan_idtones);

% group input (predictor groups)
%
grp_mtx = dummyvar(categorical(delays));
grp_mtx(:,1) = [];
grp_mtx(sum(grp_mtx,2)==0,:) = repmat(zeros(1, size(grp_mtx,2)), sum(sum(sum(grp_mtx,2)==0)), 1);
%}

%grp_mtx = delays;

% starting values
%starting_vals = [1 1 1 1 0 0 0 0];
starting_vals = [12.6 10 15 9 zeros(1, 4*size(grp_mtx,2))];


% model
[fixed_effects,re_var,stats,random_effects] = nlmefit(tone_input,wait_input,sesh_input,grp_mtx,modelFun,starting_vals);

% output     
coef_names = {'intercept', 'amplitude', 'mu', 'sigma'};
beta_pvals = 1 - tcdf( abs(fixed_effects')./stats.sebeta , stats.dfe);


% plot model
%

% colors
num_colors = size(unique(grp_mtx, 'rows'), 1);
colors = distinguishable_colors(size(unique(grp_mtx, 'rows'), 1));

% unique predictors
unq_preds = [];
for id = pos_delays
    unq_preds = [unq_preds; grp_mtx(find(delays==id,1),:)];
end

% figure on
figure; hold on
legend_trick(colors, '-')

    % plot each sessions's wait times and model
    isesh_ct = 0;
    for sesh = 1:length(delays)
        isesh_ct = isesh_ct + 1;

        % subject wait times
        if plot_on == 3
            plot(unqfrq(tone_cell{sesh}), all_wait_times{sesh}(:,1), 'o', 'color', colors(all(unq_preds == grp_mtx(isesh_ct,:),2),:))

        end

        % session model
        if ismember(plot_on, [2 3])
            PHI = fixed_effects + random_effects(:,isesh_ct);
            plot(unqfrq_interp, modelFun(PHI', all_tones, grp_mtx(isesh_ct,:)), 'color', colors(all(unq_preds == grp_mtx(isesh_ct,:),2),:))
        end

    end
    
    % plot full model
    if ismember(plot_on, [1 2 3])       
        color_ct = 0;
        for grp_predictor = 1:length(unq_preds)
            color_ct = color_ct+1;
            plot(unqfrq_interp, modelFun(fixed_effects', all_tones, unq_preds(grp_predictor,:)), '-', 'linewidth', 3, 'color', colors(color_ct,:))
        end
    
        % aesthetics
        set(gca, 'XScale', 'log')
        xlim([4600 38000])
        ylim([0 60])
        set(gca,'TickLength',[0, 0]); box off;
        rich_bounds
    end
    
    % overall plot means +/- se
    figure; hold on
    legend_trick(colors, '-')
    
    if ismember(plot_on, [1 2 3])    
        unq_delays = unique(delays);
        color_ct = 0;
        for grp_predictor = 1:length(unq_preds)
            color_ct = color_ct+1;
            
            plot(unqfrq_interp, modelFun(fixed_effects', all_tones, unq_preds(grp_predictor,:)), '-', 'color', colors(color_ct,:), 'linewidth', 4)
            plot(unqfrq_interp, modelFun(fixed_effects' + (std(random_effects(:,delays==unq_delays(grp_predictor)),[],2)./sqrt(sum(delays==unq_delays(grp_predictor))))', all_tones, unq_preds(grp_predictor,:)), '-', 'color', colors(color_ct,:), 'linewidth', 2)
            plot(unqfrq_interp, modelFun(fixed_effects' - (std(random_effects(:,delays==unq_delays(grp_predictor)),[],2)./sqrt(sum(delays==unq_delays(grp_predictor))))', all_tones, unq_preds(grp_predictor,:)), '-', 'color', colors(color_ct,:), 'linewidth', 2)
        end
        % aesthetics
        set(gca, 'XScale', 'log')
        xlim([4600 38000])
        ylim([0 60])
        set(gca,'TickLength',[0, 0]); box off;
        rich_bounds
        
        % make legend
        legend_input = []; 
        for i = 1:size(unique(grp_mtx, 'rows'), 1)
            legend_input = [legend_input; ['0' num2str(i)]];
        end
        legend(legend_input)
    end




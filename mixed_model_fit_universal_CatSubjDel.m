function [fixed_effects, random_effects, beta_pvals, coef_names, modelFun, unqfrq] = mixed_model_fit_universal_CatSubjDel(probe_paths, varargin)
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
session_ct_idx = 1:size(probe_paths,1);
%session_ct_idx = [1:10 12:13];


% Array of covariate values
    d_idx = strfind(probe_paths{1},'d');
    d_idx = d_idx(end);

    delays = [];
    for i = 1:length(probe_paths)
        
        d_idx = strfind(probe_paths{i},'d');
        d_idx = d_idx(end);
        
        probe_paths{i}(d_idx-5 : d_idx-4)
        
        delays = [delays; str2double(probe_paths{i}(d_idx-5 : d_idx-4))];
    end
    delays = delays -1 ;
    pos_delays = 0:max(delays);  

% further constrain sessions
%session_ct_idx = delays<=3;
    delays = delays(session_ct_idx);
    probe_paths = probe_paths(session_ct_idx);
    num_paths = size(probe_paths,1);
    
% print paths
probe_paths
    
    
% get unique subject ids
unq_subj = [];
for ipath = 1:size(probe_paths,1)
    gen_loc = strfind(probe_paths{ipath},'gen'); gen_loc = max(gen_loc);
    unq_subj = [unq_subj; probe_paths{ipath}(gen_loc-9 : gen_loc-2)];
end
[~, ~, unq_subj] = unique(unq_subj, 'rows', 'stable');


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
% Subject as a categorical var
%
modelFun = @(coefEsts,tones,group_mtx) sum([coefEsts(2), group_mtx.*coefEsts(6+(size(group_mtx,2)-1)*1:6+(size(group_mtx,2)-1)*2)])...
    .*exp(-((tones-sum([coefEsts(3), group_mtx.*coefEsts(7+(size(group_mtx,2)-1)*2:7+(size(group_mtx,2)-1)*3)]))...
    ./sum([coefEsts(4), group_mtx.*coefEsts(8+(size(group_mtx,2)-1)*3:8+(size(group_mtx,2)-1)*4)])).^2) + sum([coefEsts(1), group_mtx.*coefEsts(5:5+(size(group_mtx,2)-1)*1)]);
    % coefs 1:4 are {'intercept', 'amplitude', 'mu', 'sigma'}. Subsequent
    % coefs are four groups of size(group_mtx,2) coefs describing the
    % relationship between group 1 and the predictors' 'intercept's (first group), ...
    % 'amplitude's (second group), 'mu's (third group), ...
    % and 'sigma's (last group)}. 
%}




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

% group input (predictor groups)
%
%dv_del = delays;
dv_del = dummyvar(categorical(delays)); dv_del(:,1) = []
dv_subj = dummyvar(categorical(unq_subj)); dv_subj(:,1) = []
grp_mtx = [dv_del dv_subj];
%}

%grp_mtx = delays;

% starting values
%starting_vals = [20 20 20 10 zeros(1, 4*size(grp_mtx,2))];
starting_vals = [12.6 10 15 9 zeros(1, 4*size(grp_mtx,2))];


% model
[fixed_effects,re_var,stats,random_effects] = nlmefit(tone_input,wait_input,sesh_input,grp_mtx,modelFun,starting_vals);

% isolate fixed effects
fe_core = fixed_effects(1:4); % overall mean parameters
fe_rem = fixed_effects(length(fe_core)+1:end);
[fe_delay, fe_subj] = isolate_effects(fe_rem, length(fe_core), size(dv_del,2));

% internal
function [del_out, subj_out] = isolate_effects(rem_in, num_cor_params, num_del_grps)
    rem_in = reshape(rem_in,length(rem_in)/num_cor_params,num_cor_params);
    del_out = rem_in(1:num_del_grps,:);
    subj_out = rem_in(num_del_grps+1 : end,:);
end

% remove subject effects
fe_subj = zeros(size(fe_subj));

% replace fixed effects
fe_rem = [fe_delay; fe_subj];
fe_rem = fe_rem(:);
fixed_effects = [fe_core; fe_rem];


% output     
coef_names = {'intercept', 'amplitude', 'mu', 'sigma'};
beta_pvals = 1 - tcdf( abs(fixed_effects')./stats.sebeta , stats.dfe);


% plot model
%
% unique predictors
unq_preds = [];
for id = pos_delays
    unq_preds = [unq_preds; grp_mtx(find(delays==id,1),:)];
end
unq_preds(:, size(dv_del,2)+1:end) = 0;


% colors
colors = distinguishable_colors(size(unq_preds, 1));

% figure on
figure; hold on
legend_trick(colors, '-')

    % plot each sessions's wait times and model
    isesh_ct = 0;
    for sesh = 1:length(delays)
        isesh_ct = isesh_ct + 1;

        % subject wait times
        if plot_on == 3
            plot(unqfrq(tone_cell{sesh}), all_wait_times{sesh}(:,1), 'o', 'color', colors(all(unq_preds(:,1:size(dv_del,2)) == grp_mtx(isesh_ct,1:size(dv_del,2)),2),:))

        end

        % session model
        if ismember(plot_on, [2 3])
            PHI = fixed_effects + random_effects(:,isesh_ct);
            %grp_input = grp_mtx(isesh_ct,:);
            grp_input = [dv_del(isesh_ct,:) zeros(size(dv_subj(isesh_ct,:)))];
            
            %{
            unq_preds
            grp_mtx
            grp_mtx(isesh_ct,:)
            all(unq_preds(:,1:size(dv_del,2)) == grp_mtx(isesh_ct,1:size(dv_del,2)),2)
            colors(all(unq_preds(:,1:size(dv_del,2)) == grp_mtx(isesh_ct,1:size(dv_del,2)),2),:)
            %}
            
            plot(unqfrq_interp, modelFun(PHI', all_tones, grp_input), 'color', colors(all(unq_preds(:,1:size(dv_del,2)) == grp_mtx(isesh_ct,1:size(dv_del,2)),2),:))
        end

    end
    
    % plot full model
    if ismember(plot_on, [1 2 3])       

        color_ct = 0;
        for grp_predictor = 1:size(unq_preds,1)
            color_ct = color_ct+1;
            grp_input = unq_preds(grp_predictor,:);
            plot(unqfrq_interp, modelFun(fixed_effects', all_tones, grp_input), '-', 'linewidth', 3, 'color', colors(color_ct,:))
        end
    
        % aesthetics
        set(gca, 'XScale', 'log')
        xlim([4600 38000])
        ylim([0 60])
        set(gca,'TickLength',[0, 0]); box off;
        %rich_bounds
    end
    
    % overall plot means +/- se
    figure; hold on
    legend_trick(colors, '-')
    
    if ismember(plot_on, [1 2 3])    
        unq_delays = unique(delays);
        color_ct = 0;
        for grp_predictor = 1:size(unq_preds,1)
            color_ct = color_ct+1;
            %grp_input = [unq_preds(grp_predictor,:) zeros(1, size(dv_subj,2))];
            grp_input = unq_preds(grp_predictor,:);

            plot(unqfrq_interp, modelFun(fixed_effects', all_tones, grp_input), '-', 'color', colors(color_ct,:), 'linewidth', 4)
            plot(unqfrq_interp, modelFun(fixed_effects' + (std(random_effects(:,delays==unq_delays(grp_predictor)),[],2)./sqrt(sum(delays==unq_delays(grp_predictor))))', all_tones, grp_input), '-', 'color', colors(color_ct,:), 'linewidth', 2)
            plot(unqfrq_interp, modelFun(fixed_effects' - (std(random_effects(:,delays==unq_delays(grp_predictor)),[],2)./sqrt(sum(delays==unq_delays(grp_predictor))))', all_tones, grp_input), '-', 'color', colors(color_ct,:), 'linewidth', 2)
        end
        % aesthetics
        set(gca, 'XScale', 'log')
        xlim([4600 38000])
        ylim([0 60])
        set(gca,'TickLength',[0, 0]); box off;
        %rich_bounds
        
        % make legend
        legend_input = []; 
        for i = 1:size(unq_preds, 1)
            legend_input = [legend_input; ['0' num2str(i)]];
        end
        legend(legend_input)
    end

% prepare output (with all subject effects removed)
subj_idx = fixed_effects==0;
beta_pvals = beta_pvals(~subj_idx)';
random_effects = random_effects(~subj_idx,:);    
fixed_effects = fixed_effects(~subj_idx);

end

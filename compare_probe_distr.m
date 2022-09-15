function [coef_cell] = compare_probe_distr(datafolder, probe_num, day_delay, new_learn_num)
% compute normal fit statistics from probe numbers and experimental 
% conditions defined by the input 
%
% inputs are cell arrays with 1 cell per comparison group, a cell may
% contain multiple qualifiers. Empty-set inputs are ignored. All non-empty
% inputs must have the same number of cells.

% folderpath
folderpath = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\';
folderpath = [folderpath datafolder];

% check inputs
nempty_idx = [~isempty(probe_num), ~isempty(day_delay), ~isempty(new_learn_num)];
all_inp = [{probe_num},{day_delay},{new_learn_num}]; 
if any(nempty_idx)
    nempty_inp = all_inp{find(nempty_idx==1,1)};   
    num_grps = length(nempty_inp);
else
    num_grps = 0;
end
if any(cellfun('length', all_inp(nempty_idx))~=num_grps)
    error('all non-empty inputs must have same number of cells')
end


% preallocate output cells
coef_cell = cell(1, num_grps);
gof = cell(1, num_grps);
fit_p = cell(1, num_grps);
subjs = cell(1, num_grps);

%iterate through groups
figure; hold on
colors = get(gca,'ColorOrder');
for igrp = 1:num_grps
    
    % group inputs
    if ~isempty(probe_num) 
        probe_num_grp = ['-0' num2str(probe_num{igrp})]; 
    else
        probe_num_grp = [];
    end
    if ~isempty(day_delay)
        daynum_hold = num2str(day_delay{igrp});
        if length(daynum_hold)==1
            daynum_hold = ['0' daynum_hold];
        end
            day_delay_grp = [daynum_hold 'd'];
    else
        day_delay_grp = [];
    end
    if ~isempty(new_learn_num)
        new_learn_num_grp = ['learn0' num2str(new_learn_num{igrp})];
    else
        new_learn_num_grp = [];
    end
    
    % build models
    grp_paths = get_file_paths_targeted_II(folderpath, {'probe', probe_num_grp}, {day_delay_grp});
    if strcmp(day_delay_grp, '01d')
        grp_paths = grp_paths(~contains(grp_paths, '30d'));
    elseif strcmp(day_delay_grp, '30d')
        grp_paths = grp_paths(~contains(grp_paths, '01d'));
    end
    
    grp_paths
    
    
    [fixed_effects, random_effects, ~, coef_names, modelfun, unqfrq] = mixed_model_fit(grp_paths, 0);
    ses_re = std(random_effects,[],2)./sum(~isnan(random_effects),2);
    plot(unqfrq, modelfun(fixed_effects, 1:length(unqfrq)),'-','linewidth', 2, 'color', colors(igrp,:))
    plot(unqfrq, modelfun(fixed_effects+ses_re, 1:length(unqfrq)),'-','linewidth', 1, 'color', colors(igrp,:))
    plot(unqfrq, modelfun(fixed_effects-ses_re, 1:length(unqfrq)),'-','linewidth', 1, 'color', colors(igrp,:))
    set(gca, 'XScale', 'log')
    set(gca,'TickLength',[0, 0]); box off;

    % load output
    coef_cell{igrp} = fixed_effects + random_effects;
    
end

%TRY
%

figure;
for icoef = 1:length(fixed_effects)
    subplot(1,5,icoef)

    ccell = cell(1,num_grps);
    for igrp = 1:num_grps
        ccell{igrp} = coef_cell{igrp}(icoef,:)';
    end

    errorbar_plot(ccell)
    ylabel(coef_names{icoef})
end



end




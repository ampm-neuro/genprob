function [supervect, supervect_idx] = image_compare_cellType_activity(session_list, cell_regist_mtx, time_series_event_spacing)
% plot the information properties of topDown and bottomUp cells
% consider sorting inputs chronologically



%% compute trial activity for every cell in every session
%
% preallocate
session_mtx_cell = cell(size(session_list,1), 1);

% iterate through sessions
for isesh = 1:size(session_list,1)
    
    % load session
    load(session_list{isesh})

    % all trials
    hm_cell_trials = unique(trl_idx);
    
    % session matrix for every cell
    hm_cells = tw_activity_trial_hm_full_warp(trl_mtx, trl_idx, frame_times, traces, 1:size(traces,1), hm_cell_trials, [3 12], time_series_event_spacing);
    
    % remove times bins with nans on any trial
    for icell = 1:length(hm_cells)
        hm_cells{icell} = hm_cells{icell}(:,~isnan(sum(hm_cells{icell},1)));
    end

    % load
    session_mtx_cell{isesh} = hm_cells;
    
end
%save('comp_cell_type_info.mat', '-v7.3')
disp('1')
%}



%% Reorganize sessions to match universal cell registry
%
% iterate through sessions
sessions = 1:size(session_list,1);
for isesh = 1:size(session_list,1)
    
    % matrix of zeros for inactive cells
    inactive_mtx = zeros(size(session_mtx_cell{isesh}{1}));

    % session matrix with all cells
    session_mtx_cell_hold = cell(size(cell_regist_mtx,1),1);
    
    % preload all cells as inactive
    for icell = 1:length(cell_regist_mtx(:,sessions(isesh)))
        session_mtx_cell_hold{icell} = inactive_mtx;
    end
    
    % load active cells
    for icell = unique(cell_regist_mtx(cell_regist_mtx(:,sessions(isesh))>0,sessions(isesh)))'
    	session_mtx_cell_hold{find(cell_regist_mtx(:,sessions(isesh))==icell, 1)} = session_mtx_cell{isesh}{icell};
    end
    
    % load updated session matrix
    session_mtx_cell{isesh} = session_mtx_cell_hold;
    
end

%save('comp_cell_type_info.mat', '-v7.3')
disp('1.5')
%}



%% compute dispersion and reliability
%
% dispersion, not info content
[dispersion_mtx, reliability_mtx] = all_cells_infoContent(1:size(cell_regist_mtx,1), session_mtx_cell);

% set dispersion to 0 : 1
num_time_bins = 1279; %size(session_mtx_cell{1}{1},2);
divisor = 10;
dispersion_mtx = 1-dispersion_mtx./ceil(num_time_bins/divisor);

%save('comp_cell_type_info.mat', '-v7.3')
disp('2')
%}



%% identify top down and bottom up cells
%
[topDown_idx, bottomUp_idx] = cell_type_idx(cell_regist_mtx);
%save('comp_cell_type_info.mat', '-v7.3')
disp('3')
%}



%% probe and problem sessions
%{
probe_sessions = [1:3:19 20];
problem_sessions = setdiff(1:20, probe_sessions);
problem_sessions = reshape(problem_sessions, 2,6);
save('comp_cell_type_info.mat', '-v7.3')
disp('4')
%}



%% load above info
%load('comp_cell_type_info.mat')



%% top down cells
%
top_down_probes = probe_sessions(1:5);

% preallocate
tdc_disp = cell(length(top_down_probes), 3);
tdc_reli = cell(length(top_down_probes), 3);
ntdc_disp = cell(length(top_down_probes), 3);
ntdc_reli = cell(length(top_down_probes), 3);


% iterate through each probe
for iprobe = 1:length(top_down_probes)
    
    % current probe session
    cps = top_down_probes(iprobe);
    
    td_idx = topDown_idx(:,cps)>0;
    
    % probe and problem sessions
    for isesh = 1:3
        
        % current session
        cs = cps + isesh -1;
        
        % active idx
        active_idx = cell_regist_mtx(:,cs)>0;
        
        % dispersion
        tdc_disp{iprobe,isesh} = dispersion_mtx(active_idx & td_idx, cs);
        ntdc_disp{iprobe,isesh} = dispersion_mtx(active_idx & ~td_idx, cs);
        
        % reliability
        tdc_reli{iprobe,isesh} = reliability_mtx(active_idx & td_idx, cs);
        ntdc_reli{iprobe,isesh} = reliability_mtx(active_idx & ~td_idx, cs);
    end
end

figure;
for iprobe = 1:size(tdc_disp,1)
    subplot(size(tdc_disp,1),1,iprobe)
    hold on
    
    local_ebp( ntdc_disp(iprobe,:), 0.8.*[1 1 1]); % light is other
    local_ebp( tdc_disp(iprobe,:), 0.3.*[1 1 1]); % dark is coi
    
    set(gca,'TickLength',[0, 0]); box off;
    xticks(1:length(tdc_disp(iprobe,:)))
end
sgtitle('Top down dispersion')

figure;
for iprobe = 1:size(tdc_reli,1)
    subplot(size(tdc_reli,1),1,iprobe)
    hold on
    
    local_ebp( ntdc_reli(iprobe,:), 0.8.*[1 1 1]);  % light is other
    local_ebp( tdc_reli(iprobe,:), 0.3.*[1 1 1]); % dark is coi
    
    set(gca,'TickLength',[0, 0]); box off;
    xticks(1:length(tdc_reli(iprobe,:)))
end
sgtitle('Top down reliability')
%}



%% bottom up cells
%
% preallocate
buc_disp = cell(size(problem_sessions,2), 3);
buc_reli = cell(size(problem_sessions,2), 3);
nbuc_disp = cell(size(problem_sessions,2), 3);
nbuc_reli = cell(size(problem_sessions,2), 3);


% iterate through each problem
for iproblem = 1:size(problem_sessions,2)
    
    % current problems
    cur_problems = problem_sessions(:,iproblem);
    
    % bottom up index
    bu_idx = sum(bottomUp_idx(:,cur_problems),2)>0;
    
    % subsequent probe
    subs_probe = min(probe_sessions(probe_sessions>max(cur_problems)));
    
    % relevant session numbers 
    rel_sesh = [cur_problems' subs_probe];
    
    % probe and problem sessions
    for isesh = 1:length(rel_sesh)

        % current session
        cs = rel_sesh(isesh);

        % active idx
        active_idx = cell_regist_mtx(:,cs)>0;
        
        % dispersion
        buc_disp{iproblem,isesh} = dispersion_mtx(active_idx & bu_idx, cs);
        nbuc_disp{iproblem,isesh} = dispersion_mtx(active_idx & ~bu_idx, cs);
        
        % reliability
        buc_reli{iproblem,isesh} = reliability_mtx(active_idx & bu_idx, cs);
        nbuc_reli{iproblem,isesh} = reliability_mtx(active_idx & ~bu_idx, cs);
    end
end

figure;
for iproblem = 1:size(buc_disp,1)
    subplot(size(buc_disp,1),1,iproblem)
    hold on
    
    local_ebp( nbuc_disp(iproblem,:), 0.8.*[1 1 1]); % light is other
    local_ebp( buc_disp(iproblem,:), 0.3.*[1 1 1]); % dark is coi
    
    set(gca,'TickLength',[0, 0]); box off;
    xticks(1:length(buc_disp(iproblem,:)))
end
sgtitle('Bottom up dispersion')

figure;
for iproblem = 1:size(buc_reli,1)
    subplot(size(buc_reli,1),1,iproblem)
    hold on
    
    local_ebp( nbuc_reli(iproblem,:), 0.8.*[1 1 1]);  % light is other
    local_ebp( buc_reli(iproblem,:), 0.3.*[1 1 1]); % dark is coi
    
    set(gca,'TickLength',[0, 0]); box off;
    xticks(1:length(buc_reli(iproblem,:)))
end
sgtitle('Bottom up reliability')
%}



%% TOP DOWN vs New vs Other
%
% iterate through problems
for iproblem = 1:size(problem_sessions,2)

    % current problems
    cur_problems = problem_sessions(:,iproblem);
    
    % reactivated index
    reactivated_idx = sum(cell_regist_mtx(:,1:min(cur_problems)-1),2)>0;
    
    % top down
    prior_probe = max(probe_sessions(probe_sessions<min(cur_problems)));
    td_idx = topDown_idx(:,prior_probe)>0;
    
    % p8 idx
    p8_idx = cell_regist_mtx(:,20)>0;
    
    % probe and problem sessions
    for isesh = 1:length(cur_problems)

        % current session
        cs = cur_problems(isesh);

        % active idx
        active_idx = cell_regist_mtx(:,cs)>0;
        
        % dispersion
        %problemNew_disp{iproblem,isesh} = dispersion_mtx(active_idx & ~reactivated_idx, cs);
        %react_disp{iproblem,isesh} = dispersion_mtx(active_idx & td_idx, cs);
        
        % reliability
        new_reli_problem{iproblem,isesh} = (reliability_mtx(active_idx & ~reactivated_idx, cs));
        topDown_reli_problem{iproblem,isesh} = (reliability_mtx(active_idx & td_idx, cs));
        otherRet_reli_problem{iproblem,isesh} = (reliability_mtx(active_idx & ~td_idx & reactivated_idx, cs));

    end
end

% RELIABILITY behavior learning curve style fig
figure; hold on
for iprob = 1:size(new_reli_problem,1)
    local_ebp(problem_sessions(:,iprob),  new_reli_problem(iprob,:), 0.8.*[1 1 1]);  % light is new
    local_ebp(problem_sessions(:,iprob),  topDown_reli_problem(iprob,:), 0.3.*[1 1 1]); % dark is coi
    local_ebp(problem_sessions(:,iprob),  otherRet_reli_problem(iprob,:), [.1 .4 .8]); % color is other
end
set(gca,'TickLength',[0, 0]); box off;
xticks(2.5:3:20)
xticklabels(1:6)
xlabel('Problem')
xlabel('Reliability')
title('Reliability in new and retrieved neurons')

% RELIABILITY combining across problems fig
%[a b c d] = ttest2(cell2mat(new_reli(:,1)), cell2mat(topDown_reli(:,1)))
new_reli_problem = [{cell2mat(new_reli_problem(:,1))} {cell2mat(new_reli_problem(:,2))}];
topDown_reli_problem = [{cell2mat(topDown_reli_problem(:,1))} {cell2mat(topDown_reli_problem(:,2))}];
otherRet_reli_problem = [{cell2mat(otherRet_reli_problem(:,1))} {cell2mat(otherRet_reli_problem(:,2))}];
figure; hold on
    local_ebp([1 4], new_reli_problem, 0.8.*[1 1 1]);  % light is new
    local_ebp([1.75 4.75], topDown_reli_problem, 0.3.*[1 1 1]); % dark is coi
    local_ebp([2.5 5.5], otherRet_reli_problem, [.1 .4 .8]); % color is other
set(gca,'TickLength',[0, 0]); box off;
xticks(1:2)
xticklabels({'First', 'Last'})
xlabel('Problem')
ylabel('Reliability')
title('Reliability in new and retrieved neurons')
set(gcf, 'Position', [29   596   424   745])
%}



%% BOTTOM UP vs New vs Other

% iterate through problems
for iproblem = 1:size(problem_sessions,2)
    
    % current problems
    cur_problems = problem_sessions(:,iproblem);
    
    % bottom up
    next_probe = min(probe_sessions(probe_sessions>max(cur_problems)));
    bu_idx = sum(bottomUp_idx(:,cur_problems),2)>0;
    
    % reactivated index
    reactivated_idx = sum(cell_regist_mtx(:,1:next_probe-1),2)>0;
    
    % active idx
    active_idx = cell_regist_mtx(:,next_probe)>0;

    % reliability
    new_reli_probe{iproblem} = reliability_mtx(active_idx & ~reactivated_idx, next_probe);
    bottomUp_reli_probe{iproblem} = reliability_mtx(active_idx & bu_idx, next_probe);
    otherRet_reli_probe{iproblem} = reliability_mtx(active_idx & ~bu_idx & reactivated_idx, next_probe);
    
end
new_reli_probe = new_reli_probe';
bottomUp_reli_probe = bottomUp_reli_probe';
otherRet_reli_probe = otherRet_reli_probe';

% RELIABILITY behavior learning curve style fig
figure; hold on
for iprob = 1:size(new_reli_probe,1)
    local_ebp(probe_sessions(iprob),  new_reli_probe(iprob,:), 0.8.*[1 1 1]);  % light is new
    local_ebp(probe_sessions(iprob)+0.75,  bottomUp_reli_probe(iprob,:), 0.3.*[1 1 1]); % dark is coi
    local_ebp(probe_sessions(iprob)+1.50,  otherRet_reli_probe(iprob,:), [.1 .4 .8]); % color is other
end
set(gca,'TickLength',[0, 0]); box off;
%xticks(2.5:3:20)
%xticklabels(1:6)
xlabel('Probe')
ylabel('Reliability')
title('Reliability in new and retrieved neurons')

% RELIABILITY combining across problems fig
new_reli_probe = {cell2mat(new_reli_probe)};
bottomUp_reli_probe = {cell2mat(bottomUp_reli_probe)};
otherRet_reli_probe = {cell2mat(otherRet_reli_probe)};
figure; hold on
    local_ebp(1, new_reli_probe, 0.8.*[1 1 1]);  % light is new
    local_ebp(2, bottomUp_reli_probe, 0.3.*[1 1 1]); % dark is coi
    local_ebp(3, otherRet_reli_probe, [.1 .4 .8]); % color is other
set(gca,'TickLength',[0, 0]); box off;
xticks(1:2)
xticklabels({'New', 'Bottom up', 'Other reactivated'})
ylabel('Reliability')
title('Reliability in new and retrieved neurons')



% bar plot combining problem and probe cells
supervect = [new_reli_problem{1}; topDown_reli_problem{1}; otherRet_reli_problem{1}; ...
    new_reli_problem{2}; topDown_reli_problem{2}; otherRet_reli_problem{2}; ...
    new_reli_probe{:}; bottomUp_reli_probe; otherRet_reli_probe]';
supervect = cell2mat(supervect');
supervect_idx = [zeros(size(new_reli_problem{1})); ones(size(topDown_reli_problem{1})); 2.*ones(size(otherRet_reli_problem{1})); ...
    3.*ones(size(new_reli_problem{2})); 4.*ones(size(topDown_reli_problem{2})); 5.*ones(size(otherRet_reli_problem{2})); ...
    6.*ones(size(cell2mat(new_reli_probe))); 7.*ones(size(cell2mat(bottomUp_reli_probe))); 8.*ones(size(cell2mat(otherRet_reli_probe)))];
figure; hold on
%boxplot(supervect, supervect_idx)
%
    local_ebp([1.00 4.00], new_reli_problem, 0.8.*[1 1 1]);  % light is new
    local_ebp([1.75 4.75], topDown_reli_problem, 0.3.*[1 1 1]); % dark is coi
    local_ebp([2.50 5.50], otherRet_reli_problem, [.1 .4 .8]); % color is other
    
    local_ebp(7.00, new_reli_probe, 0.8.*[1 1 1]);  % light is new
    local_ebp(7.75, bottomUp_reli_probe, 0.3.*[1 1 1]); % dark is coi
    local_ebp(8.50, otherRet_reli_probe, [.1 .4 .8]); % color is other
    set(gca,'TickLength',[0, 0]); box off;
    ylabel('Reliability')
    xlabel('Session')
    xticks(1.75:3:10)
    xticklabels({'First problem','Last problem', 'Probe'})
%}


    
end

function local_ebp(varargin)


if length(varargin)==3
    xpos = varargin{1};
    cell_e = varargin{2};
    plot_color = varargin{3};
elseif length(varargin)==2
    cell_e = varargin{1};
    plot_color = varargin{2};
    xpos = 1:length(cell_e);
end

means = nan(length(cell_e),1);
stds = nan(length(cell_e),1);
sqrt_l = nan(length(cell_e),1);
cell_e_xaxis = cell(size(cell_e));

for imtx = 1: length(cell_e)
    
        mean_imtx = nanmean(cell_e{imtx});
        std_imtx = max([nanstd(cell_e{imtx}) realmin]);
        for idpt = 1:length(cell_e{imtx})
           
            base_jitter_multiplier = length(cell_e{imtx})*0.0006 + 0.375;
            base_jitter = (rand(1)-0.5)*base_jitter_multiplier;
            dist_from_mean = abs(cell_e{imtx}(idpt) - mean_imtx);
            std_from_mean = dist_from_mean/std_imtx;
            bulb_correction = std_from_mean/7.5 + (rand(1)-0.5)*0.4;
            cell_e_xaxis{imtx}(idpt) = xpos(imtx)+base_jitter*(1-bulb_correction);
            plot(cell_e_xaxis{imtx}(idpt), cell_e{imtx}(idpt), 'o', 'color', plot_color, 'markersize', 5)
        end
        
    means(imtx) = nanmean(cell_e{imtx});
    stds(imtx) = nanstd(cell_e{imtx});
    sqrt_l(imtx) = sqrt(sum(~isnan(cell_e{imtx})));

end
std_es = stds./sqrt_l;

errorbar(xpos, means, std_es, 'k.', 'linewidth', 1.5)
plot(xpos, means, 'color', plot_color, 'linewidth', 1.5)
%bar(xpos, means)



end
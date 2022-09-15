function [abr_ON,abp_ON,abr_OFF,abp_OFF,asubj,aprobs] = probe_change_opto(probes_in)
% evaluates how the laser ON and laser OFF responses on the current probe
% compare to the laser OFF responses on the previous probe

% opto path
opto_root_path = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_optoExp_acc';

% rich and poor tone indices
load('unqfrq41', 'unqfrq41')
all_prob_tones = rich_bounds_prob('mevar', 0);
for iprob = 1:size(all_prob_tones,1)
   tones_hold = find(ismember(unqfrq41,all_prob_tones(iprob,:)));
   all_prob_tones(iprob,1) = tones_hold(ismember(tones_hold, 16:26));
   all_prob_tones(iprob,2) = tones_hold(ismember(tones_hold, setdiff(1:41, 16:26))); 
end

% preallocate
unq_subj_cell = cell(1,length(probes_in));
means_ON = cell(size(unq_subj_cell));
std_ON = cell(size(unq_subj_cell));
se_ON = cell(size(unq_subj_cell));
means_OFF = cell(size(unq_subj_cell));
std_OFF = cell(size(unq_subj_cell));
se_OFF = cell(size(unq_subj_cell));

% get all probe responses
for iprobe = 1:length(probes_in)
    
    % get paths
    if probes_in(iprobe) <=6
        path_str_1 = 'preprobe_opto';
        path_str_2 = ['_0' num2str(probes_in(iprobe)) '_'];
    elseif probes_in(iprobe) > 6
        path_str_1 = 'postprobe_opto';
        path_str_2 = ['_0' num2str(probes_in(iprobe)-6) '_'];
    end
    paths = get_file_paths_targeted(opto_root_path, path_str_1, path_str_2);
    
    % load unique subjects
    subj_list = [];
    for ipath = 1:size(paths,1)
        subj_list = [subj_list; find_subj_id(paths{ipath})];
    end
    unq_subj_cell{iprobe} = unique(subj_list, 'rows');
    
    % get session responses
    [mean_wait_times_optoON, std_wait_times_optoON, se_wait_times_optoON,...
        mean_wait_times_optoOFF, std_wait_times_optoOFF, se_wait_times_optoOFF]...
        = plot_allprobes_opto_smooth(paths, 0);

    % load responses
    means_ON{iprobe} = mean_wait_times_optoON;
    std_ON{iprobe} = std_wait_times_optoON;
    se_ON{iprobe} = se_wait_times_optoON;
    means_OFF{iprobe} = mean_wait_times_optoOFF;
    std_OFF{iprobe} = std_wait_times_optoOFF;
    se_OFF{iprobe} = se_wait_times_optoOFF;
    
end

% preallocate
all_change_ON = cell(size(unq_subj_cell));
all_change_OFF = cell(size(unq_subj_cell));
beyond_rich_change_ON = cell(size(unq_subj_cell));
beyond_poor_change_ON = cell(size(unq_subj_cell));
beyond_rich_change_OFF = cell(size(unq_subj_cell));
beyond_poor_change_OFF = cell(size(unq_subj_cell));
all_shared_subjs = cell(size(unq_subj_cell));

% compute probe changes
for iprobe = 1:length(probes_in)-1
    
    % shared subjects indices
    last_subjs = unq_subj_cell{iprobe};
    cur_subjs = unq_subj_cell{iprobe+1};
    [all_shared_subjs{iprobe}, last_subj_idx, cur_subj_idx] = ...
        intersect(last_subjs, cur_subjs, 'rows', 'stable');
    
    % last probe
    last_means = means_OFF{iprobe}(last_subj_idx,:);
    
    % current probe
    cur_means_ON = means_ON{iprobe+1}(cur_subj_idx,:);
    cur_means_OFF = means_OFF{iprobe+1}(cur_subj_idx,:);
    
    % changes
    all_change_ON{iprobe} = (cur_means_ON - last_means)./(cur_means_ON + last_means);
    all_change_OFF{iprobe} = (cur_means_OFF - last_means)./(cur_means_OFF + last_means);
    
    % beyond rich and poor changes
    rich_tone = all_prob_tones(probes_in(iprobe), 1);
    poor_tone = all_prob_tones(probes_in(iprobe), 2);
    if rich_tone > poor_tone
        beyond_rich_change_ON{iprobe} = all_change_ON{iprobe}(:,rich_tone:end);
        beyond_poor_change_ON{iprobe} = all_change_ON{iprobe}(:,1:poor_tone);
        beyond_rich_change_OFF{iprobe} = all_change_OFF{iprobe}(:,rich_tone:end);
        beyond_poor_change_OFF{iprobe} = all_change_OFF{iprobe}(:,1:poor_tone);
    else
        beyond_rich_change_ON{iprobe} = all_change_ON{iprobe}(:,1:rich_tone);
        beyond_poor_change_ON{iprobe} = all_change_ON{iprobe}(:,poor_tone:end);
        beyond_rich_change_OFF{iprobe} = all_change_OFF{iprobe}(:,1:rich_tone);
        beyond_poor_change_OFF{iprobe} = all_change_OFF{iprobe}(:,poor_tone:end);
    end
    
    % plot smoothed lines
    %
    figure; hold on
    
        % ON
        plot(mean(all_change_ON{iprobe}), 'b-', 'linewidth', 2)
        plot(mean(all_change_ON{iprobe})-std(all_change_ON{iprobe})./sqrt(size(all_change_ON{iprobe},1)), 'b-')
        plot(mean(all_change_ON{iprobe})+std(all_change_ON{iprobe})./sqrt(size(all_change_ON{iprobe},1)), 'b-')

        % OFF
        plot(mean(all_change_OFF{iprobe}), 'k-', 'linewidth', 2)
        plot(mean(all_change_OFF{iprobe})-std(all_change_OFF{iprobe})./sqrt(size(all_change_OFF{iprobe},1)), 'k-')
        plot(mean(all_change_OFF{iprobe})+std(all_change_OFF{iprobe})./sqrt(size(all_change_OFF{iprobe},1)), 'k-')
    
    plot(xlim, [0 0], 'k--')
    set(gca,'TickLength',[0, 0]); box off;
    ylim([-.4 .4])
    title(['Change from probe ' num2str(iprobe) ' to probe ' num2str(iprobe+1)])
    ylabel('Change in wait times (normalized)')
    xlabel('Frequency (Khz)')
    frq_disp = 1:5:length(unqfrq41);
    xticks(frq_disp); xticklabels(unqfrq41(frq_disp)./1000)
    xlim([0 42])
        
        % problem tones
        plot(rich_tone.*[1 1], ylim, 'r-')
        plot(poor_tone.*[1 1], ylim, '-', 'color', 0.4.*[1 1 1])
    
end

% plot
abr_ON = [];
abp_ON = [];
abr_OFF = [];
abp_OFF = [];
asubj = [];
aprobs = [];
for ib = 1:length(beyond_rich_change_ON)
    abr_ON = [abr_ON; mean(beyond_rich_change_ON{ib},2)];
    abp_ON = [abp_ON; mean(beyond_poor_change_ON{ib},2)];
    abr_OFF = [abr_OFF; mean(beyond_rich_change_OFF{ib},2)];
    abp_OFF = [abp_OFF; mean(beyond_poor_change_OFF{ib},2)];
    asubj = [asubj; all_shared_subjs{ib}];
    aprobs = [aprobs; repmat(ib, size(mean(beyond_rich_change_ON{ib},2)))];
end
figure; errorbar_barplot([{abp_OFF}, {abp_ON}, {abr_OFF}, {abr_ON}], 1);
xticks(1:4)
xticklabels({'poor;OFF', 'poor;ON', 'rich;OFF', 'rich;ON'})
title('Probe-over-probe change')
ylabel('Wait times')

% model beyond changes
%
change_dv = [abr_ON abp_ON abr_OFF abp_OFF];
opto_idx = [ones(size(abr_ON)) ones(size(abp_ON)) zeros(size(abr_OFF)) zeros(size(abp_OFF))];
rich_idx = [ones(size(abr_ON)) zeros(size(abp_ON)) ones(size(abr_OFF)) zeros(size(abp_OFF))];
    [~,~,unq_subj_nums] = unique(asubj, 'rows'); 
subject_idx = repmat(unq_subj_nums, 1, size(change_dv,2));
problem_idx = repmat(aprobs, 1, size(change_dv,2));

tbl = table(change_dv(:), opto_idx(:), rich_idx(:), subject_idx(:), problem_idx(:),...
    'VariableNames',{'change','opto','rich','subj','problem'});
tbl.opto = categorical(tbl.opto);
tbl.rich = categorical(tbl.rich);
tbl.subject = categorical(tbl.subj);
tbl.problem = categorical(tbl.problem);

model_str = 'change~opto+rich+opto*rich+(1|subject)+(1|problem)';
lme = fitlme(tbl, model_str);

[~,p_poor,~,s_poor] = ttest(abp_OFF, abp_ON)
[~,p_rich,~,stats_rich] = ttest(abr_OFF, abr_ON)






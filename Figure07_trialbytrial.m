green_color = [21 154 72]./255; blue_color =  [43 73 159]./255;
light_green_color = [162 201 129]./255; light_blue_color = [102 164 213]./255;
figure_path = 'C:\Users\ampm1\Documents\manuscripts\gentrain\manuscript_figures\fig06\';
beh_fld = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\';
cell_reg_path = 'C:\Users\ampm1\Documents\MATLAB\generalization\cellreg_data\';
predict_subjs = [{'651049m1'} {'658648m2'} {'691359m1'} {'690330m1'} {'691359m2'} {'696944m3'} {'699437m4'} {'699438m3'}];
unpredict_subjs = [{'687034m4'} {'683472m2'} {'699437m2'} {'699438m1'} {'683472m3'}];
session_chron_reorder_idx = [1 9 15 2 10 16 3 11 17 4 12 18 5 13 19 6 14 20 7 8];

% for time warping
tses = [0.2 1.0 1.0 2.0 3.0 2.0];
event_frame = cumsum(tses(1:end-1)).*100; %from second np to reward delivery/ quit
event_frame_rng = event_frame(1):event_frame(4); % NP2 to end of min wait

% key imaging sessions
preprobe_session_numbers = 1:3:16;
firstday_session_numbers = 2:3:17;
lastday_session_numbers = 3:3:18;

% SEE ALSO: t_by_t_learning


%% compute firing activity for every neuron on every first-session trial

% predict
%{
subj_cell_activity_firstsesh_predict = cell(1, length(predict_subjs));
for isubj = 1:length(predict_subjs)
    subj_cell_activity_firstsesh_predict{isubj} = cell(1,6);
    for iproblem = 1:6
         load(problem_paths_predict{isubj}{iproblem})
         subj_cell_activity_firstsesh_predict{isubj}{iproblem} = tw_activity_plot_trial_hm_warp(trl_mtx, trl_idx, frame_times, traces, 1:size(traces,1), 1:size(trl_mtx,1), ceil([sum(tses(1:3)) sum(tses(4:end))]), tses);
         subj_cell_activity_firstsesh_predict{isubj}{iproblem} = cat(3, subj_cell_activity_firstsesh_predict{isubj}{iproblem}{:}); %trials, time, neuron
    end
end
save('first_session_activity_predict', 'subj_cell_activity_firstsesh_predict')

% unpredict
subj_cell_activity_firstsesh_unpredict = cell(1, length(unpredict_subjs));
for isubj = 1:length(unpredict_subjs)
    subj_cell_activity_firstsesh_unpredict{isubj} = cell(1,6);
    for iproblem = 1:6
         load(problem_paths_unpredict{isubj}{iproblem})
         subj_cell_activity_firstsesh_unpredict{isubj}{iproblem} = tw_activity_plot_trial_hm_warp(trl_mtx, trl_idx, frame_times, traces, 1:size(traces,1), 1:size(trl_mtx,1), ceil([sum(tses(1:3)) sum(tses(4:end))]), tses);
         subj_cell_activity_firstsesh_unpredict{isubj}{iproblem} = cat(3, subj_cell_activity_firstsesh_unpredict{isubj}{iproblem}{:}); %trials, time, neuron
    end
end
save('first_session_activity_unpredict', 'subj_cell_activity_firstsesh_unpredict')
%}
load('first_session_activity_predict', 'subj_cell_activity_firstsesh_predict')
load('first_session_activity_unpredict', 'subj_cell_activity_firstsesh_unpredict')



%% constrain to desired event frame range
%{
% predict
subj_cell_activity_firstsesh_predict_constrained = subj_cell_activity_firstsesh_predict;
for isubj = 1:length(predict_subjs)
    for iproblem = 1:6
        all_nan_idx = all(isnan(subj_cell_activity_firstsesh_predict_constrained{isubj}{iproblem}(:,:,1)));
        subj_cell_activity_firstsesh_predict_constrained{isubj}{iproblem} = subj_cell_activity_firstsesh_predict_constrained{isubj}{iproblem}(:,~all_nan_idx,:);
        subj_cell_activity_firstsesh_predict_constrained{isubj}{iproblem} = subj_cell_activity_firstsesh_predict_constrained{isubj}{iproblem}(:,event_frame_rng,:);
    end
end
% unpredict
subj_cell_activity_firstsesh_unpredict_constrained = subj_cell_activity_firstsesh_unpredict;
for isubj = 1:length(unpredict_subjs)
    for iproblem = 1:6
        all_nan_idx = all(isnan(subj_cell_activity_firstsesh_unpredict_constrained{isubj}{iproblem}(:,:,1)));
        subj_cell_activity_firstsesh_unpredict_constrained{isubj}{iproblem} = subj_cell_activity_firstsesh_unpredict_constrained{isubj}{iproblem}(:,~all_nan_idx,:);
        subj_cell_activity_firstsesh_unpredict_constrained{isubj}{iproblem} = subj_cell_activity_firstsesh_unpredict_constrained{isubj}{iproblem}(:,event_frame_rng,:);
    end
end
%save
save('first_session_activity_constrained', 'subj_cell_activity_firstsesh_predict_constrained', 'subj_cell_activity_firstsesh_unpredict_constrained')
%}
load('first_session_activity_constrained', 'subj_cell_activity_firstsesh_predict_constrained', 'subj_cell_activity_firstsesh_unpredict_constrained')



%% compute trial by trial (adjacent) change in correlation
%{ 
% predict
adjacent_corrs_predict = cell(length(predict_subjs), 6);
for isubj = 1:length(subj_cell_activity_firstsesh_predict_constrained)
    for iproblem = 1:length(subj_cell_activity_firstsesh_predict_constrained{isubj})
        [isubj iproblem]
        adjacent_corrs_predict{isubj, iproblem} = ...
            nan(size(subj_cell_activity_firstsesh_predict_constrained{isubj}{iproblem},3), size(subj_cell_activity_firstsesh_predict_constrained{isubj}{iproblem},1)); %neuron, trial comp
        for ic = 1:size(subj_cell_activity_firstsesh_predict_constrained{isubj}{iproblem},3)
            for itrl = 1:size(subj_cell_activity_firstsesh_predict_constrained{isubj}{iproblem},1)-1
                adjacent_corrs_predict{isubj, iproblem}(ic, itrl) =...
                    corr(subj_cell_activity_firstsesh_predict_constrained{isubj}{iproblem}(itrl,:,ic)', subj_cell_activity_firstsesh_predict_constrained{isubj}{iproblem}(itrl+1,:,ic)');
            end
        end
    end
end

% unpredict
adjacent_corrs_unpredict = cell(length(unpredict_subjs), 6);
for isubj = 1:length(subj_cell_activity_firstsesh_unpredict_constrained)
    for iproblem = 1:length(subj_cell_activity_firstsesh_unpredict_constrained{isubj})
        [isubj iproblem]
        adjacent_corrs_unpredict{isubj, iproblem} = ...
            nan(size(subj_cell_activity_firstsesh_unpredict_constrained{isubj}{iproblem},3), size(subj_cell_activity_firstsesh_unpredict_constrained{isubj}{iproblem},1)); %neuron, trial comp
        for ic = 1:size(subj_cell_activity_firstsesh_unpredict_constrained{isubj}{iproblem},3)
            for itrl = 1:size(subj_cell_activity_firstsesh_unpredict_constrained{isubj}{iproblem},1)-1
                adjacent_corrs_unpredict{isubj, iproblem}(ic, itrl) =...
                    corr(subj_cell_activity_firstsesh_unpredict_constrained{isubj}{iproblem}(itrl,:,ic)', subj_cell_activity_firstsesh_unpredict_constrained{isubj}{iproblem}(itrl+1,:,ic)');
            end
        end
    end
end
save('first_session_adjacent_trial_corrs', 'adjacent_corrs_predict', 'adjacent_corrs_unpredict')
%}
load('first_session_adjacent_trial_corrs', 'adjacent_corrs_predict', 'adjacent_corrs_unpredict')



%% compute mean of adjacent corr matrices
%{
% predict
mean_acm_predict = cell(length(predict_subjs), 6);
for isubj = 1:size(adjacent_corrs_predict,1)
    for iproblem = 1:size(adjacent_corrs_predict,2)
        mean_acm_predict{isubj, iproblem} = mean(adjacent_corrs_predict{isubj, iproblem});
    end
end

% unpredict
mean_acm_unpredict = cell(length(unpredict_subjs), 6);
for isubj = 1:size(adjacent_corrs_unpredict,1)
    for iproblem = 1:size(adjacent_corrs_unpredict,2)
        mean_acm_unpredict{isubj, iproblem} = mean(adjacent_corrs_unpredict{isubj, iproblem});
    end
end
%}

%% interp and plot adjacent means

%{
% predict
mean_acm_predict_interp = nan(length(predict_subjs), 600);
for isubj = 1:size(adjacent_corrs_predict,1)
    iprob_rng = 1:100;
    for iproblem = 1:size(adjacent_corrs_predict,2)
        raw_vals = mean_acm_predict{isubj, iproblem};
        mean_acm_predict_interp(isubj, iprob_rng) = interp1(1:length(raw_vals), raw_vals, linspace(1, length(raw_vals), 100));
        iprob_rng = iprob_rng+100;
    end
end

% unpredict
mean_acm_unpredict_interp = nan(length(unpredict_subjs), 100);
for isubj = 1:size(adjacent_corrs_unpredict,1)
    iprob_rng = 1:100;
    for iproblem = 1:size(adjacent_corrs_unpredict,2)
        raw_vals = mean_acm_unpredict{isubj, iproblem};
        mean_acm_unpredict_interp(isubj, iprob_rng) = interp1(1:length(raw_vals), raw_vals, linspace(1, length(raw_vals), 100));
        iprob_rng = iprob_rng+100;
    end
end

% interped figure
figure; hold on

    % unpredict
    for isubj = 1:size(mean_acm_unpredict_interp,1)
        %plot(mean_acm_unpredict_interp(isubj,:), 'color', light_blue_color, 'linewidth', 1)
    end
    % predict
    for isubj = 1:size(mean_acm_predict_interp,1)
        %plot(mean_acm_predict_interp(isubj,:), 'color', light_green_color, 'linewidth', 1)
    end
    
    % unpredict
    all_corr_mean = nanmean(mean_acm_unpredict_interp);
    all_corr_ste = nanstd(mean_acm_unpredict_interp)./sqrt(sum(~isnan(mean_acm_unpredict_interp(:,1))));
    plot(all_corr_mean, 'color', blue_color, 'linewidth', 2.5)
    plot(all_corr_mean-all_corr_ste, 'color', blue_color, 'linewidth', 1)
    plot(all_corr_mean+all_corr_ste, 'color', blue_color, 'linewidth', 1)

    % predict
    all_corr_mean = nanmean(mean_acm_predict_interp);
    all_corr_ste = nanstd(mean_acm_predict_interp)./sqrt(sum(~isnan(mean_acm_predict_interp(:,1))));
    plot(all_corr_mean, 'color', green_color, 'linewidth', 2.5)
    plot(all_corr_mean-all_corr_ste, 'color', green_color, 'linewidth', 1)
    plot(all_corr_mean+all_corr_ste, 'color', green_color, 'linewidth', 1)
    
    % problem bounds
    ylim([-.6 .6]); set(gca,'TickLength',[0, 0]);
    plot([100.5:100:600]'.*[1 1], ylim, 'k-', 'linewidth', 2)
    % no preference
    plot(xlim, [0 0], 'k--', 'linewidth', 2)
    
    % test fig
    figure; fit_line([suprise_mtx_predict(:); suprise_mtx_unpredict(:)], [mean_acm_predict_interp(:); mean_acm_unpredict_interp(:)])
    xlabel('Instantaneous Suprise')
    ylabel('Instantaneous Stability')
%}

drawnow
%% compute mean pre-probe activity

% probe paths
%{
    % predict
    preprobe_paths_predict = cell(length(predict_subjs), 6);
    for isubj = 1:length(predict_subjs)
        preprobe_paths_predict{isubj} = get_file_paths_targeted(beh_fld, {predict_subjs{isubj}, 'preprobe'});
    end
    
    % unpredict
    preprobe_paths_unpredict = cell(length(unpredict_subjs), 6);
    for isubj = 1:length(unpredict_subjs)
        preprobe_paths_unpredict{isubj} = get_file_paths_targeted(beh_fld, {unpredict_subjs{isubj}, 'preprobe'});
    end

    
% compute activity
%
    % predict
    subj_cell_activity_preprobe_predict = cell(1, length(predict_subjs));
    for isubj = 1:length(predict_subjs)
        subj_cell_activity_preprobe_predict{isubj} = cell(1,6);
        for iproblem = 1:6
            [isubj iproblem]
             load(preprobe_paths_predict{isubj}{iproblem})
             subj_cell_activity_preprobe_predict{isubj}{iproblem} = tw_activity_plot_trial_hm_warp(trl_mtx, trl_idx, frame_times, traces, 1:size(traces,1), 1:size(trl_mtx,1), ceil([sum(tses(1:3)) sum(tses(4:end))]), tses);
             subj_cell_activity_preprobe_predict{isubj}{iproblem} = cat(3, subj_cell_activity_preprobe_predict{isubj}{iproblem}{:}); %trials, time, neuron
        end
    end
    save('preprobe_session_activity_predict', 'subj_cell_activity_preprobe_predict')

    % unpredict
    subj_cell_activity_preprobe_unpredict = cell(1, length(unpredict_subjs));
    for isubj = 1:length(unpredict_subjs)
        subj_cell_activity_preprobe_unpredict{isubj} = cell(1,6);
        for iproblem = 1:6
            [isubj iproblem]
             load(preprobe_paths_unpredict{isubj}{iproblem})
             subj_cell_activity_preprobe_unpredict{isubj}{iproblem} = tw_activity_plot_trial_hm_warp(trl_mtx, trl_idx, frame_times, traces, 1:size(traces,1), 1:size(trl_mtx,1), ceil([sum(tses(1:3)) sum(tses(4:end))]), tses);
             subj_cell_activity_preprobe_unpredict{isubj}{iproblem} = cat(3, subj_cell_activity_preprobe_unpredict{isubj}{iproblem}{:}); %trials, time, neuron
        end
    end
    save('preprobe_session_activity_unpredict', 'subj_cell_activity_preprobe_unpredict')
    %}
load('preprobe_session_activity_predict', 'subj_cell_activity_preprobe_predict')
load('preprobe_session_activity_unpredict', 'subj_cell_activity_preprobe_unpredict')



%% average and constrain preprobe sessions to desired event frame range
%{
% predict
subj_cell_activity_preprobe_predict_constrained = subj_cell_activity_preprobe_predict;
for isubj = 1:length(predict_subjs)
    for iprobe = 1:6
        all_nan_idx = all(isnan(subj_cell_activity_preprobe_predict_constrained{isubj}{iprobe}(:,:,1)));
        subj_cell_activity_preprobe_predict_constrained{isubj}{iprobe} = subj_cell_activity_preprobe_predict_constrained{isubj}{iprobe}(:,~all_nan_idx,:);
        subj_cell_activity_preprobe_predict_constrained{isubj}{iprobe} = subj_cell_activity_preprobe_predict_constrained{isubj}{iprobe}(:,event_frame_rng,:);
        subj_cell_activity_preprobe_predict_constrained{isubj}{iprobe} = permute(mean(subj_cell_activity_preprobe_predict_constrained{isubj}{iprobe}), [3 2 1]);
    end
end
% unpredict
subj_cell_activity_preprobe_unpredict_constrained = subj_cell_activity_preprobe_unpredict;
for isubj = 1:length(unpredict_subjs)
    for iprobe = 1:6
        all_nan_idx = all(isnan(subj_cell_activity_preprobe_unpredict_constrained{isubj}{iprobe}(:,:,1)));
        subj_cell_activity_preprobe_unpredict_constrained{isubj}{iprobe} = subj_cell_activity_preprobe_unpredict_constrained{isubj}{iprobe}(:,~all_nan_idx,:);
        subj_cell_activity_preprobe_unpredict_constrained{isubj}{iprobe} = subj_cell_activity_preprobe_unpredict_constrained{isubj}{iprobe}(:,event_frame_rng,:);
        subj_cell_activity_preprobe_unpredict_constrained{isubj}{iprobe} = permute(mean(subj_cell_activity_preprobe_unpredict_constrained{isubj}{iprobe}), [3 2 1]);
    end
end
%save
save('preprobe_activity_constrained', 'subj_cell_activity_preprobe_predict_constrained', 'subj_cell_activity_preprobe_unpredict_constrained')
%}
load('preprobe_activity_constrained', 'subj_cell_activity_preprobe_predict_constrained', 'subj_cell_activity_preprobe_unpredict_constrained')



%% isolate held cells using CellReg output
%{
% predict
subj_cell_activity_preprobe_predict_constrained_held = subj_cell_activity_preprobe_predict_constrained;
subj_cell_activity_firstsesh_predict_constrained_held = subj_cell_activity_firstsesh_predict_constrained;
for isubj = 1:length(predict_subjs)

    % load cell reg
    load([cell_reg_path predict_subjs{isubj} '\cell_reg_' predict_subjs{isubj} '.mat'], 'cell_registered_struct')   
    if strcmp(predict_subjs{isubj}, '690330m1')
        cell_reg_idx = cell_registered_struct.cell_to_index_map;
        cell_reg_idx = [cell_reg_idx(:,1:9) nan(size(cell_reg_idx(:,1))) cell_reg_idx(:,10:end)];
    else
        cell_reg_idx = cell_registered_struct.cell_to_index_map;
    end
    cell_reg_idx = cell_reg_idx(:,session_chron_reorder_idx);

    % iterate through probe-problem pairs
    for iprob = 1:6
    
        % held cell index
        preprobe_firstday_nums = [preprobe_session_numbers(iprob) firstday_session_numbers(iprob)];
        held_cell_idx_probe = cell_reg_idx(all(cell_reg_idx(:,preprobe_firstday_nums)>0,2), preprobe_firstday_nums(1));
        held_cell_idx_problem = cell_reg_idx(all(cell_reg_idx(:,preprobe_firstday_nums)>0,2), preprobe_firstday_nums(2));
    
        % only keep held cells
        %
            % probe
            subj_cell_activity_preprobe_predict_constrained_held{isubj}{iprob} = ...
                subj_cell_activity_preprobe_predict_constrained_held{isubj}{iprob}(held_cell_idx_probe, :);
            % problem
            subj_cell_activity_firstsesh_predict_constrained_held{isubj}{iprob} = ...
                subj_cell_activity_firstsesh_predict_constrained_held{isubj}{iprob}(:, :, held_cell_idx_problem);
    end
end

% unpredict
subj_cell_activity_preprobe_unpredict_constrained_held = subj_cell_activity_preprobe_unpredict_constrained;
subj_cell_activity_firstsesh_unpredict_constrained_held = subj_cell_activity_firstsesh_unpredict_constrained;
for isubj = 1:length(unpredict_subjs)

    % load cell reg
    load([cell_reg_path unpredict_subjs{isubj} '\cell_reg_' unpredict_subjs{isubj} '.mat'], 'cell_registered_struct')
    cell_reg_idx = cell_registered_struct.cell_to_index_map;
    cell_reg_idx = cell_reg_idx(:,session_chron_reorder_idx);
    
    % iterate through probe-problem pairs
    for iprob = 1:6
    
        % held cell index
        preprobe_firstday_nums = [preprobe_session_numbers(iprob) firstday_session_numbers(iprob)];
        held_cell_idx_probe = cell_reg_idx(all(cell_reg_idx(:,preprobe_firstday_nums)>0,2), preprobe_firstday_nums(1));
        held_cell_idx_problem = cell_reg_idx(all(cell_reg_idx(:,preprobe_firstday_nums)>0,2), preprobe_firstday_nums(2));
    
        % only keep held cells
        %
            % probe
            subj_cell_activity_preprobe_unpredict_constrained_held{isubj}{iprob} = ...
                subj_cell_activity_preprobe_unpredict_constrained_held{isubj}{iprob}(held_cell_idx_probe, :);
            % problem
            subj_cell_activity_firstsesh_unpredict_constrained_held{isubj}{iprob} = ...
                subj_cell_activity_firstsesh_unpredict_constrained_held{isubj}{iprob}(:, :, held_cell_idx_problem);
    end
end
save('held_cell_activity', 'subj_cell_activity_preprobe_unpredict_constrained_held', 'subj_cell_activity_firstsesh_unpredict_constrained_held')
%}
load('held_cell_activity', 'subj_cell_activity_preprobe_unpredict_constrained_held', 'subj_cell_activity_firstsesh_unpredict_constrained_held')



%% compute problem trial by problem trial correlation to preprobe for every held cell
%{
% predict
probe_corrs_predict = cell(length(predict_subjs), 6);
for isubj = 1:length(subj_cell_activity_firstsesh_predict_constrained_held)
    for iprob = 1:length(subj_cell_activity_firstsesh_predict_constrained_held{isubj})
        probe_corrs_predict{isubj, iprob} = ...
            nan(size(subj_cell_activity_firstsesh_predict_constrained_held{isubj}{iprob},3), size(subj_cell_activity_firstsesh_predict_constrained_held{isubj}{iprob},1)); %neuron, trial comp
        for ic = 1:size(subj_cell_activity_firstsesh_predict_constrained_held{isubj}{iprob},3)
            for itrl = 1:size(subj_cell_activity_firstsesh_predict_constrained_held{isubj}{iprob},1)
                probe_mean_activity = subj_cell_activity_preprobe_predict_constrained_held{isubj}{iprob}(ic,:)';
                problem_trial_activity = subj_cell_activity_firstsesh_predict_constrained_held{isubj}{iprob}(itrl,:,ic)';
                probe_corrs_predict{isubj, iprob}(ic, itrl) = corr(probe_mean_activity, problem_trial_activity);
            end
        end
    end
end

% unpredict
probe_corrs_unpredict = cell(length(unpredict_subjs), 6);
for isubj = 1:length(subj_cell_activity_firstsesh_unpredict_constrained_held)
    for iprob = 1:length(subj_cell_activity_firstsesh_unpredict_constrained_held{isubj})
        probe_corrs_unpredict{isubj, iprob} = ...
            nan(size(subj_cell_activity_firstsesh_unpredict_constrained_held{isubj}{iprob},3), size(subj_cell_activity_firstsesh_unpredict_constrained_held{isubj}{iprob},1)); %neuron, trial comp
        for ic = 1:size(subj_cell_activity_firstsesh_unpredict_constrained_held{isubj}{iprob},3)
            for itrl = 1:size(subj_cell_activity_firstsesh_unpredict_constrained_held{isubj}{iprob},1)
                probe_mean_activity = subj_cell_activity_preprobe_unpredict_constrained_held{isubj}{iprob}(ic,:)';
                problem_trial_activity = subj_cell_activity_firstsesh_unpredict_constrained_held{isubj}{iprob}(itrl,:,ic)';
                probe_corrs_unpredict{isubj, iprob}(ic, itrl) = corr(probe_mean_activity, problem_trial_activity);
            end
        end
    end
end
save('preprobe_corrs', 'probe_corrs_predict', 'probe_corrs_unpredict')
%}
load('preprobe_corrs', 'probe_corrs_predict', 'probe_corrs_unpredict')

%% compute mean of probe corr matrices
%

% minimum number of cells for valid mean correlation (try 5)
min_cells = 5;
probe_corrs_cellct_predict = nan(length(predict_subjs), 6);
probe_corrs_cellct_unpredict = nan(length(unpredict_subjs), 6);

% predict
mean_acm_probe_predict = cell(length(predict_subjs), 6);
for isubj = 1:size(probe_corrs_predict,1)
    for iproblem = 1:size(probe_corrs_predict,2)
        probe_corrs_cellct_predict(isubj, iproblem) = size(probe_corrs_predict{isubj, iproblem},1);
        if probe_corrs_cellct_predict(isubj, iproblem)>=min_cells
            mean_acm_probe_predict{isubj, iproblem} = mean(probe_corrs_predict{isubj, iproblem},1);
        else
            mean_acm_probe_predict{isubj, iproblem} = nan(size(mean(probe_corrs_predict{isubj, iproblem},1)));
        end
    end
end

% unpredict
mean_acm_probe_unpredict = cell(length(unpredict_subjs), 6);
for isubj = 1:size(probe_corrs_unpredict,1)
    for iproblem = 1:size(probe_corrs_unpredict,2)
        probe_corrs_cellct_unpredict(isubj, iproblem) = size(probe_corrs_unpredict{isubj, iproblem},1);
        if probe_corrs_cellct_unpredict(isubj, iproblem)>=min_cells
            mean_acm_probe_unpredict{isubj, iproblem} = mean(probe_corrs_unpredict{isubj, iproblem},1);
        else
            mean_acm_probe_unpredict{isubj, iproblem} = nan(size(mean(probe_corrs_unpredict{isubj, iproblem},1)));
        end
    end
end
%}



%% interp and plot probe-to-problem correlation means

%{
% predict
mean_acm_probe_predict_interp = nan(length(predict_subjs), 600);
for isubj = 1:size(probe_corrs_predict,1)
    iprob_rng = 1:100;
    for iproblem = 1:size(probe_corrs_predict,2)
        raw_vals = mean_acm_probe_predict{isubj, iproblem};
        %if all_raw_vals
        mean_acm_probe_predict_interp(isubj, iprob_rng) = interp1(1:length(raw_vals), raw_vals, linspace(1, length(raw_vals), 100));
        iprob_rng = iprob_rng+100;
    end
end

% unpredict
mean_acm_probe_unpredict_interp = nan(length(unpredict_subjs), 100);
for isubj = 1:size(probe_corrs_unpredict,1)
    iprob_rng = 1:100;
    for iproblem = 1:size(probe_corrs_unpredict,2)
        raw_vals = mean_acm_probe_unpredict{isubj, iproblem};
        mean_acm_probe_unpredict_interp(isubj, iprob_rng) = interp1(1:length(raw_vals), raw_vals, linspace(1, length(raw_vals), 100));
        iprob_rng = iprob_rng+100;
    end
end
save('Stability_correlations_interped', 'mean_acm_probe_predict_interp', 'mean_acm_probe_unpredict_interp')
%}
load('Stability_correlations_interped', 'mean_acm_probe_predict_interp', 'mean_acm_probe_unpredict_interp')

% interped figure
figure; hold on

    % unpredict
    for isubj = 1:size(mean_acm_probe_unpredict_interp,1)
        %plot(mean_acm_probe_unpredict_interp(isubj,:), 'color', light_blue_color, 'linewidth', 1)
    end
    % predict
    for isubj = 1:size(mean_acm_probe_predict_interp,1)
        %plot(mean_acm_probe_predict_interp(isubj,:), 'color', light_green_color, 'linewidth', 1)
    end
    
    % unpredict
    all_corr_mean = nanmean(mean_acm_probe_unpredict_interp);
    all_corr_ste = nanstd(mean_acm_probe_unpredict_interp)./sqrt(sum(~isnan(mean_acm_probe_unpredict_interp(:,1))));
    plot(all_corr_mean, 'color', blue_color, 'linewidth', 2.5)
    plot(all_corr_mean-all_corr_ste, 'color', blue_color, 'linewidth', 1)
    plot(all_corr_mean+all_corr_ste, 'color', blue_color, 'linewidth', 1)

    % predict
    all_corr_mean = nanmean(mean_acm_probe_predict_interp);
    all_corr_ste = nanstd(mean_acm_probe_predict_interp)./sqrt(sum(~isnan(mean_acm_probe_predict_interp(:,1))));
    plot(all_corr_mean, 'color', green_color, 'linewidth', 2.5)
    plot(all_corr_mean-all_corr_ste, 'color', green_color, 'linewidth', 1)
    plot(all_corr_mean+all_corr_ste, 'color', green_color, 'linewidth', 1)
    
    % problem bounds
    ylim([-1 1]); set(gca,'TickLength',[0, 0]);
    plot([100.5:100:600]'.*[1 1], ylim, 'k-', 'linewidth', 2)
    % no preference
    plot(xlim, [0 0], 'k--', 'linewidth', 2)
    title( 'stability')

    
    
%% compute suprise scores for every trial
%
[suprise_cell_predict, discrim_cell_predict, suprise_mtx_predict, problem_paths_predict] = suprise(predict_subjs');
[suprise_cell_unpredict, discrim_cell_unpredict, suprise_mtx_unpredict, problem_paths_unpredict] = suprise(unpredict_subjs');
%save('suprise_scores', 'suprise_cell_predict', 'discrim_cell_predict','suprise_mtx_predict', 'problem_paths_predict',...
%    'suprise_cell_unpredict', 'discrim_cell_unpredict', 'suprise_mtx_unpredict', 'problem_paths_unpredict');
%}
%load('suprise_scores', 'suprise_cell_predict', 'discrim_cell_predict','suprise_mtx_predict', 'problem_paths_predict',...
%    'suprise_cell_unpredict', 'discrim_cell_unpredict', 'suprise_mtx_unpredict', 'problem_paths_unpredict');


% remove sessions with too few neurons for estimating mean correlations
%{
problem_range = 1:100;
for iproblem = 1:size(probe_corrs_cellct_predict,2)
    suprise_mtx_predict(probe_corrs_cellct_predict(:,iproblem)<min_cells, problem_range) = nan;
    suprise_mtx_unpredict(probe_corrs_cellct_unpredict(:,iproblem)<min_cells, problem_range) = nan;
    problem_range = problem_range+100;
end
%}

% interped figure
figure; hold on

    % unpredict
    for isubj = 1:size(suprise_mtx_unpredict,1)
        plot(suprise_mtx_unpredict(isubj,:), 'color', light_blue_color, 'linewidth', 1)
    end
    % predict
    for isubj = 1:size(suprise_mtx_predict,1)
        plot(suprise_mtx_predict(isubj,:), 'color', light_green_color, 'linewidth', 1)
    end
    
    % unpredict
    all_supr_mean = nanmean(suprise_mtx_unpredict);
    all_supr_ste = nanstd(suprise_mtx_unpredict)./sqrt(sum(~isnan(suprise_mtx_unpredict)));
    plot(all_supr_mean, 'color', blue_color, 'linewidth', 2.5)
    plot(all_supr_mean-all_supr_ste, 'color', blue_color, 'linewidth', 1)
    plot(all_supr_mean+all_supr_ste, 'color', blue_color, 'linewidth', 1)

    % predict
    all_supr_mean = nanmean(suprise_mtx_predict);
    all_supr_ste = nanstd(suprise_mtx_predict)./sqrt(sum(~isnan(suprise_mtx_predict)));
    plot(all_supr_mean, 'color', green_color, 'linewidth', 2.5)
    plot(all_supr_mean-all_supr_ste, 'color', green_color, 'linewidth', 1)
    plot(all_supr_mean+all_supr_ste, 'color', green_color, 'linewidth', 1)
    
    % problem bounds
    ylim([-1 1]); set(gca,'TickLength',[0, 0]);
    plot([100.5:100:600]'.*[1 1], ylim, 'k-', 'linewidth', 2)
    % no preference
    plot(xlim, [0 0], 'k--', 'linewidth', 2)
    title Suprise
    
        
    
%% psuedo session

% psuedo

    % suprise
        figure; hold on
        errorbar_mtx(nanmean((reshape(suprise_mtx_unpredict,[5, 100, 6])),3))
        errorbar_mtx(nanmean((reshape(suprise_mtx_predict,[8, 100, 6])),3))
        plot(xlim, [0 0], 'k--');
    % stability
        figure; hold on
        errorbar_mtx(nanmean(reshape((mean_acm_probe_unpredict_interp')',[5, 100, 6]),3))
        errorbar_mtx(nanmean(reshape((mean_acm_probe_predict_interp')',[8, 100, 6]),3))
        plot(xlim, [0 0], 'k--');
        
        
%% trial pace

% predict
%{
subj_cell_initiationtime_firstsesh_predict = cell(1, length(predict_subjs));
for isubj = 1:length(predict_subjs)
    subj_cell_initiationtime_firstsesh_predict{isubj} = cell(1,6);
    for iproblem = 1:6
         load(problem_paths_predict{isubj}{iproblem})
         subj_cell_initiationtime_firstsesh_predict{isubj}{iproblem} = abs(trl_mtx(:,6)); %NP on to end of fixed delay (length trials)
    end
end

% unpredict
subj_cell_initiationtime_firstsesh_unpredict = cell(1, length(unpredict_subjs));
for isubj = 1:length(unpredict_subjs)
    subj_cell_initiationtime_firstsesh_unpredict{isubj} = cell(1,6);
    for iproblem = 1:6
         load(problem_paths_unpredict{isubj}{iproblem})
         subj_cell_initiationtime_firstsesh_unpredict{isubj}{iproblem} = abs(trl_mtx(:,6)); %NP on to end of fixed delay (length trials)
    end
end
save('trial_initation_times', 'subj_cell_initiationtime_firstsesh_predict', 'subj_cell_initiationtime_firstsesh_unpredict')
%}
load('trial_initation_times', 'subj_cell_initiationtime_firstsesh_predict', 'subj_cell_initiationtime_firstsesh_unpredict')


%% trial start time
% predict
%{
subj_cell_starttime_firstsesh_predict = cell(1, length(predict_subjs));
for isubj = 1:length(predict_subjs)
    subj_cell_starttime_firstsesh_predict{isubj} = cell(1,6);
    for iproblem = 1:6
         load(problem_paths_predict{isubj}{iproblem})
         subj_cell_starttime_firstsesh_predict{isubj}{iproblem} = abs(trl_mtx(:,1)); %NP on to end of fixed delay (length trials)
    end
end

% unpredict
subj_cell_starttime_firstsesh_unpredict = cell(1, length(unpredict_subjs));
for isubj = 1:length(unpredict_subjs)
    subj_cell_starttime_firstsesh_unpredict{isubj} = cell(1,6);
    for iproblem = 1:6
         load(problem_paths_unpredict{isubj}{iproblem})
         subj_cell_starttime_firstsesh_unpredict{isubj}{iproblem} = abs(trl_mtx(:,1)); %NP on to end of fixed delay (length trials)
    end
end
save('trial_start_times', 'subj_cell_starttime_firstsesh_predict', 'subj_cell_starttime_firstsesh_unpredict')
%}
load('trial_start_times', 'subj_cell_starttime_firstsesh_predict', 'subj_cell_starttime_firstsesh_unpredict')




%% plot 2 var figure

% colormap
green_colors = [162 221 154; 118 205 117; 065 181 093; 037 149 069; 007 120 056; 013 078 032]./255;
blue_colors = [158 202 225; 106 173 214; 067 146 198; 035 114 181; 009 083 157; 033 052 104]./255;

% open figure
    figure; hold on

    
    % dot plot 
    prob_cells_unpredict = cell(2,6);
    prob_cells_predict = cell(2,6);
    for iprobe = 1:length(suprise_cell_unpredict{1})
        %unpredict
        for isubj = 1:length(suprise_cell_unpredict)
                plot(suprise_cell_unpredict{isubj}{iprobe}, mean_acm_probe_unpredict{isubj, iprobe}', 'o', 'color', blue_colors(iprobe,:))
                prob_cells_unpredict{1,iprobe} = [prob_cells_unpredict{1,iprobe}; suprise_cell_unpredict{isubj}{iprobe}]; % suprise
                prob_cells_unpredict{2,iprobe} = [prob_cells_unpredict{2,iprobe}; mean_acm_probe_unpredict{isubj, iprobe}']; % stability
        end
        
        % predict
        for isubj = 1:length(suprise_cell_predict)
                plot(suprise_cell_predict{isubj}{iprobe}, mean_acm_probe_predict{isubj, iprobe}', 'o', 'color', green_colors(iprobe,:))
                prob_cells_predict{1,iprobe} = [prob_cells_predict{1,iprobe};suprise_cell_predict{isubj}{iprobe}];
                prob_cells_predict{2,iprobe} = [prob_cells_predict{2,iprobe};mean_acm_probe_predict{isubj, iprobe}'];
        end
    end    
    % means 
        % unpredict
        for iprobe = 1:length(suprise_cell_unpredict{1})
            plot(nanmean(prob_cells_unpredict{1,iprobe}), nanmean(prob_cells_unpredict{2,iprobe}), '.', 'color', blue_colors(iprobe,:), 'markersize', 60)
        % predict
            plot(nanmean(prob_cells_predict{1,iprobe}), nanmean(prob_cells_predict{2,iprobe}), '.', 'color', green_colors(iprobe,:), 'markersize', 60)
        end
    
    [r, p] = fit_line([cell2mat(prob_cells_predict(1,:)'); cell2mat(prob_cells_unpredict(1,:)')], [cell2mat(prob_cells_predict(2,:)'); cell2mat(prob_cells_unpredict(2,:)')], 0);
    title(['Stability and suprise' '; r=' num2str(r) ', p=' num2str(p)])
    xlabel('Surprise')
    ylabel('Stability')
    set(gca,'TickLength',[0, 0]); box off;
    axis square; axis([-1 1 -1 1])    
    plot(xlim, [0 0], 'k--');plot([0 0], ylim, 'k--');
    
    
    
    
%% model

% extra info
%
    group_cell_predict = cell(1,6); group_cell_unpredict = cell(1,6);
    subject_cell_predict = cell(1,6); subject_cell_unpredict = cell(1,6);
    session_cell_predict = cell(1,6); session_cell_unpredict = cell(1,6);
    session_trial_number_predict = cell(1,6); session_trial_number_unpredict = cell(1,6);
    initiationtime_predict = cell(1,6); initiationtime_unpredict = cell(1,6);
    start_time_predict = cell(1,6); start_time_unpredict = cell(1,6);
    engagement_predict = cell(1,6); engagement_unpredict = cell(1,6);
    for iprobe = 1:length(suprise_cell_unpredict{1})
        
        % predict
        for isubj_predict = 1:length(suprise_cell_predict)
            local_size = size(suprise_cell_predict{isubj_predict}{iprobe});
            group_cell_predict{1,iprobe} = [group_cell_predict{1,iprobe}; zeros(local_size)];
            subject_cell_predict{1,iprobe} = [subject_cell_predict{1,iprobe}; repmat(isubj_predict, local_size)];
            session_cell_predict{1,iprobe} = [session_cell_predict{1,iprobe}; repmat(iprobe, local_size)];
            session_trial_number_predict{1,iprobe} = [session_trial_number_predict{1,iprobe}; (1:local_size(1))'];
            initiationtime_predict{1,iprobe} = [initiationtime_predict{1,iprobe}; subj_cell_initiationtime_firstsesh_predict{isubj_predict}{iprobe}];
            start_time_predict{1,iprobe} = [start_time_predict{1,iprobe}; subj_cell_starttime_firstsesh_predict{isubj_predict}{iprobe}];
            engagement_predict{1,iprobe} = [engagement_predict{1,iprobe}; subj_cell_starttime_firstsesh_predict{isubj_predict}{iprobe}(1:end)-[0;subj_cell_starttime_firstsesh_predict{isubj_predict}{iprobe}(1:end-1)]];
        end
        
        %unpredict
        for isubj_unpredict = 1:length(suprise_cell_unpredict)
            local_size = size(suprise_cell_unpredict{isubj_unpredict}{iprobe});
            group_cell_unpredict{1,iprobe} = [group_cell_unpredict{1,iprobe}; zeros(local_size)];
            subject_cell_unpredict{1,iprobe} = [subject_cell_unpredict{1,iprobe}; repmat(isubj_predict+isubj_unpredict, local_size)];
            session_cell_unpredict{1,iprobe} = [session_cell_unpredict{1,iprobe}; repmat(iprobe, local_size)];
            session_trial_number_unpredict{1,iprobe} = [session_trial_number_unpredict{1,iprobe}; (1:local_size(1))'];
            initiationtime_unpredict{1,iprobe} = [initiationtime_unpredict{1,iprobe}; subj_cell_initiationtime_firstsesh_unpredict{isubj_unpredict}{iprobe}];
            start_time_unpredict{1,iprobe} = [start_time_unpredict{1,iprobe}; subj_cell_starttime_firstsesh_unpredict{isubj_unpredict}{iprobe}];
            engagement_unpredict{1,iprobe} = [engagement_unpredict{1,iprobe}; subj_cell_starttime_firstsesh_unpredict{isubj_unpredict}{iprobe}(1:end)-[0;subj_cell_starttime_firstsesh_unpredict{isubj_unpredict}{iprobe}(1:end-1)]];
        end
    end

stability_vect = [cell2mat(prob_cells_predict(2,:)'); cell2mat(prob_cells_unpredict(2,:)')];
suprise_vect = [cell2mat(prob_cells_predict(1,:)'); cell2mat(prob_cells_unpredict(1,:)')];
group_vect = [cell2mat(group_cell_predict'); cell2mat(group_cell_unpredict')];
session_vect = [cell2mat(session_cell_predict'); cell2mat(session_cell_unpredict')];
session_trial_num_vect = [cell2mat(session_trial_number_predict'); cell2mat(session_trial_number_unpredict')];
subject_id_vect = [cell2mat(subject_cell_predict'); cell2mat(subject_cell_unpredict')];
trial_initiation_time_vect = [cell2mat(initiationtime_predict'); cell2mat(initiationtime_unpredict')];
start_time_vect = [cell2mat(start_time_predict'); cell2mat(start_time_unpredict')];
engagment_vect = [cell2mat(engagement_predict'); cell2mat(engagement_unpredict')];

% mixed model 
%{
model_str = 'Stability~Suprise+Trial+Pace+(1|Subject)+(1|Session)';
tbl = table(stability_vect, suprise_vect, session_trial_num_vect, trial_initiation_time_vect, subject_id_vect, session_vect, 'VariableNames', {'Stability', 'Suprise', 'Trial', 'Pace', 'Subject', 'Session'});
tbl.Subject = categorical(tbl.Subject);
lme_PreprobeFinalprobe = fitlme(tbl, model_str)
%}

model_str = 'Stability~CumSuprise+NegPace+NegEngagement+(1|Subject)';
tbl = table(stability_vect, suprise_vect, trial_initiation_time_vect, start_time_vect, engagment_vect, session_trial_num_vect, group_vect, session_vect, subject_id_vect, 'VariableNames', {'Stability', 'CumSuprise', 'NegPace', 'StartTime', 'NegEngagement', 'TrialNum', 'Group', 'StageNum', 'Subject'});
tbl.Subject = categorical(tbl.Subject);
lme_PreprobeFinalprobe = fitlme(tbl, model_str)






function [first_sesh_rates_sorted, second_sesh_rates_sorted, cell_corrs, subj_idx] = pairwise_session_remapping(mevar_or_hivar, subject_ids, session_numbers)
% plots firing rate mapping of cells in the earlier session and the later
% session, using every pairwise comparision of the sessions in
% session_numbers (1 through 20)

% pairwise_session_remapping('mevar', {'651049m1','658648m2'}, [1:4])
% pairwise_session_remapping('hivar', {'683472m2','683472m3'}, [1:4])

% see pairwise_session_remapping_wrap

%% details
%tses = [0.2 1.1 1.0 2.0 8.5 2.0];
tses = [0.3 1.0 2.0 2.0 2.0 2.0];
%tses = [2 2 2 2 2 2];
%tses = [1 3 3 .5 .5 2];
stage_bins = cumsum(tses.*100);
session_chron_reorder_idx = [1 9 15 2 10 16 3 11 17 4 12 18 5 13 19 6 14 20 7 8];



%% pairwise comparisons

session_numbers = sort(session_numbers);
pairwise_comps = [];
for isesh1 = 1:length(session_numbers)
    for isesh2 = 1:length(session_numbers)
        if isesh1<isesh2
            pairwise_comps = [pairwise_comps; session_numbers([isesh1 isesh2])];
        end
    end
end



%% get session paths
subj_sesh_cell = cell(1,length(subject_ids));
for isubj = 1:length(subject_ids)
    
    % last problem session paths...
    last_prob_sessions = [];
    for iprob = 1:6
        prob_sessions = get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_' mevar_or_hivar '_imaging_hpc\' subject_ids{isubj} ], ['var0' num2str(iprob)], 'LED');
        last_prob_sessions = [last_prob_sessions; prob_sessions(end)];
    end
    
    % ...merge with pre and postprobe and first problem session paths
    session_cell = [...
        get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_' mevar_or_hivar '_imaging_hpc\' subject_ids{isubj} ], 'preprobe', 'LED');...
        get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_' mevar_or_hivar '_imaging_hpc\' subject_ids{isubj} ] , 'postprobe', 'LED');...
        get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_' mevar_or_hivar '_imaging_hpc\' subject_ids{isubj} ], 'var0', '-01', 'LED');...
        last_prob_sessions];
    
    % sort chronologically
    subj_sesh_cell{isubj} = session_cell(session_chron_reorder_idx);
end



%% load neuron IDs
subj_crm = cell(1,length(subject_ids));
for isubj = 1:length(subject_ids)
    
    % cell_regist_mtx
    clearvars cell_regist_mtx
    load(['C:\Users\ampm1\Documents\MATLAB\generalization\cellreg_data\' subject_ids{isubj} '\cell_reg_' subject_ids{isubj} '.mat']); 
    cell_regist_mtx = cell_registered_struct.cell_to_index_map;
    
    % correct for missing sessions
    if strcmp(subject_ids{isubj}, '690330m1')
        cell_regist_mtx = [cell_regist_mtx(:,1:9) nan(size(cell_regist_mtx(:,1))) cell_regist_mtx(:,10:end)];
    end
    
    % chron
    cell_regist_mtx = cell_regist_mtx(:, session_chron_reorder_idx);
    
    % load sessions of interest only
    subj_crm{isubj} = cell_regist_mtx;
    
end



%% compute mean firing rate maps for each session
mfr_maps = cell(length(session_numbers), length(subject_ids));
for isubj = 1:length(subject_ids)
    for isesh = 1:length(session_numbers)
    
        % load session
        subj_sesh_cell{isubj}{session_numbers(isesh)}
        load(subj_sesh_cell{isubj}{session_numbers(isesh)}, 'trl_mtx', 'trl_idx', 'frame_times', 'traces')

        % compute firing rates
        [~, ~, ~, ~, all_mtx_unsort] = image_mean_activity_timewarp(trl_mtx, trl_idx, frame_times, traces, 1:size(traces,1), unique(trl_idx), tses);

        
        %figure; imagesc(all_mtx_unsort); title one; drawnow; still good
        
        % only until end of random delay
        %all_mtx_unsort = all_mtx_unsort(:,1:stage_bins(4));
        all_mtx_unsort = norm_mtx(all_mtx_unsort')';
        
        % load
        mfr_maps{isesh, isubj} = all_mtx_unsort;

    end
end
    
    
    
%% load mean firing rate maps
first_sesh_rates = [];
second_sesh_rates = [];
subj_idx = [];
for isubj = 1:length(subject_ids)
    for icomp = 1:size(pairwise_comps,1)
        
        % common cells
        cell_regist_mtx_local = subj_crm{isubj}(:, pairwise_comps(icomp,:));
        cell_regist_mtx_local = cell_regist_mtx_local(sum(cell_regist_mtx_local>0,2)==2,:);
        
        sum(cell_regist_mtx_local)
        
        for isesh = 1:2
            
            
            % common cells only
            all_mtx_unsort = mfr_maps{session_numbers==pairwise_comps(icomp,isesh), isubj}(cell_regist_mtx_local(:,isesh), :);
            
            % load rates
            if isesh == 1
                first_sesh_rates = [first_sesh_rates; all_mtx_unsort];
            elseif isesh == 2
                second_sesh_rates = [second_sesh_rates; all_mtx_unsort];
            end
            subj_idx = [subj_idx; repmat(isubj, size(all_mtx_unsort,1), 1)];
            
        end
    end
end


%figure; imagesc(first_sesh_rates); title two; drawnow
%figure; imagesc(second_sesh_rates); title three; drawnow


%% sort
[first_sesh_rates_sorted, sort_idx] = sort_rows_by_peak(first_sesh_rates);
second_sesh_rates_sorted = second_sesh_rates(sort_idx,:);


%% compute correlations
cell_corrs = nan(size(first_sesh_rates_sorted,1),1);
for ic = 1:size(first_sesh_rates_sorted,1)
    cell_corrs(ic) = corr(first_sesh_rates_sorted(ic,:)', second_sesh_rates_sorted(ic,:)');
end

%% plot

% plot corms
figure;

subplot(1,2,1);
imagesc(first_sesh_rates_sorted)
set(gca,'TickLength',[0, 0]); box off;
xticks([])
hold on
for irl = 1:5
    plot(stage_bins(irl).*[1 1], ylim, 'r-')
end
title('First session')

subplot(1,2,2);
imagesc(second_sesh_rates_sorted)
set(gca,'TickLength',[0, 0]); box off;
xticks([])
hold on
for irl = 1:5
    plot(stage_bins(irl).*[1 1], ylim, 'r-')
end
title('Second session')

sgtitle([mevar_or_hivar ' ' num2str(session_numbers)])

% plot errorbar
figure; hold on
colors = distinguishable_colors(length(subject_ids));
errorbar_plot_lineonly({atanh(cell_corrs)})
tic_vect = [-.99 -.96 -.9 -.8 -.6 -.3 0 .3 .6 .8 .9 .96 .99 .999]; 
ylim(atanh([-.6 .9])); yticks(atanh(tic_vect)); yticklabels(tic_vect); ylim(atanh([-.95 .995]))
plot(xlim, [0 0], 'k--')
title([mevar_or_hivar ' ' num2str(session_numbers)])










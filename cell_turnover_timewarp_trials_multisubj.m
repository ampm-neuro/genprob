function [all_cell_corrs, all_merge_mtx, subj_corr_cell, all_common_cell_matrices, subj_corr_means_mtx, subj_corr_cell_cellid] = cell_turnover_timewarp_trials_multisubj(mevar_or_hivar, subject_ids, sessions, tses)
% runs cell_turnover_timewarp_trials on each subject and then combines
% outputs into multisubj plots

% mevar or hivar is a string 'mevar' or 'hivar'
% subject ids is a cell vector is 1 subj id per cell eg [{'651049m1'} {'658648m2'}]

% for time warping
session_chron_reorder_idx = [1 9 15 2 10 16 3 11 17 4 12 18 5 13 19 6 14 20 7 8];



%% iterate through subjects
subj_corr_cell = cell(1,length(subject_ids));
subj_corr_cell_cellid = cell(size(subj_corr_cell));
all_cell_corrs = cell(20);
all_common_cell_matrices = cell(19,2);
all_merge_mtx = [];
subj_corr_means_mtx = [];
for isubj = 1:length(subject_ids)
    subject_ids{isubj}
    
    % load cell registration matrix
    load(['cellreg_data\' subject_ids{isubj} '\cell_reg_' subject_ids{isubj} '.mat'])
    cell_regist_mtx = cell_registered_struct.cell_to_index_map;
    
    % chron
    if strcmp(subject_ids{isubj}, '690330m1')
        cell_regist_mtx = [cell_regist_mtx(:,1:9) nan(size(cell_regist_mtx(:,1))) cell_regist_mtx(:,10:end)]; 
        scri = session_chron_reorder_idx;
    elseif size(cell_regist_mtx,2)==19
        scri = session_chron_reorder_idx(1:end-1); scri(scri>=8) = scri(scri>=8)-1;
    else
        scri = session_chron_reorder_idx;
    end
    
    % desired sessions only
    scri = scri(sessions);
    cell_regist_mtx = cell_regist_mtx(:, scri);
    
    % get session paths
    last_prob_sessions = [];
    for iprob = 1:6
        prob_sessions = get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_' mevar_or_hivar '_imaging_hpc\' subject_ids{isubj} ], ['var0' num2str(iprob)], 'LED');
        last_prob_sessions = [last_prob_sessions; prob_sessions(end)];
    end
    session_cell = [...
        get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_' mevar_or_hivar '_imaging_hpc\' subject_ids{isubj} ], 'preprobe', 'LED');...
        get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_' mevar_or_hivar '_imaging_hpc\' subject_ids{isubj} ] , 'postprobe', 'LED');...
        get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_' mevar_or_hivar '_imaging_hpc\' subject_ids{isubj} ], 'var0', '-01', 'LED');...
        last_prob_sessions];
    session_cell = session_cell(scri);

    % constrain sessions
    %session_cell = session_cell(1:4);
    %cell_regist_mtx = cell_regist_mtx(:, 1:4);
    %session_cell = session_cell(17:end);
    %cell_regist_mtx = cell_regist_mtx(:, 17:end);
    
    % compute common cell correlations
    [common_cell_corrs_cell, common_cell_corrs_cell_means, merge_mtx, common_cell_matrices, common_cell_corrs_cell_cellids] = cell_turnover_timewarp_trials(session_cell, cell_regist_mtx, tses);
    subj_corr_cell{isubj} = common_cell_corrs_cell;
    subj_corr_cell_cellid{isubj} = common_cell_corrs_cell_cellids;
    subj_corr_means_mtx = cat(3, subj_corr_means_mtx, common_cell_corrs_cell_means);

    close all
    
    %{
    if strcmp(subject_ids{isubj}, '651049m1')
        common_cell_corrs_cell_temp = cell(20);
        common_cell_corrs_cell_temp(16:20, 16:20) = common_cell_corrs_cell;
        common_cell_corrs_cell = common_cell_corrs_cell_temp;
    end
    %}
    
    % combine correlations
    for i1=1:size(common_cell_corrs_cell,1)
        for i2=1:size(common_cell_corrs_cell,2)
            all_cell_corrs{i1,i2} = [all_cell_corrs{i1,i2}; common_cell_corrs_cell{i1,i2}];
            
            if i2 == i1+1
                all_common_cell_matrices{i1, 1} = [all_common_cell_matrices{i1, 1}; common_cell_matrices{i1, 1}];
                all_common_cell_matrices{i1, 2} = [all_common_cell_matrices{i1, 2}; common_cell_matrices{i1, 2}];
            end
        end
    end
    
    % combine merged matrices
    all_merge_mtx = [all_merge_mtx; merge_mtx];

end

%save('all_cell_corrs.mat', '-v7.3')


%% plots


%correlation between sessions using mutually active cells

% preallocate
common_cell_corrs_means = nan(20);

% iterate through session comparisons
for i1 = 1:size(all_cell_corrs,1)
    for i2 = 1:size(all_cell_corrs,2)
        common_cell_corrs_means(i1, i2) = nanmean(all_cell_corrs{i1, i2});
    end
end

% plot means
figure;
ampm_pcolor(common_cell_corrs_means)
axis square
caxis([-1 1])
colorbar
xlabel('Probe'); ylabel('Probe')
title(remove_underscore('cell_corrs_mutual_means'))



% errorbar line plot of probe over probe similarity

% just probe over probe
cell_corrs_mutual_temp = all_cell_corrs(2:end,1:end-1); 
cell_corrs_mutual_temp = cell_corrs_mutual_temp(logical(eye(size(cell_corrs_mutual_temp,1))));

% transform
for icell = 1:length(cell_corrs_mutual_temp)
    cell_corrs_mutual_temp{icell} = atanh(cell_corrs_mutual_temp{icell});
end

% plot
figure; errorbar_plot(cell_corrs_mutual_temp);

% aesthetics
hold on; plot(xlim, [0 0], 'k--')
ylim([-1 1])
maxr = .99; tic_vect = [-maxr -.9 -.8 -.6 -.3 0 .3 .6 .8 .9 maxr];
ylim(atanh([-maxr maxr])); 
ylim_hold = ylim;
ylim([ylim_hold(1) ylim_hold(2)+(ylim_hold(2)*.25)])
yticks(atanh(tic_vect)); yticklabels(tic_vect)


%% sorting heat maps of adjacent sessions

zoom_portion = 20:430; %end np to end fixed delay

% iterate through comparisons
for icomp = 1:size(all_common_cell_matrices,1)
    
    
    if size(all_common_cell_matrices{icomp,1},1)<1
        continue
    end
    
    % sort first session
    [~, sort_idx] = sort_rows_by_peak(all_common_cell_matrices{icomp,1}(:,zoom_portion));
    
    % plot first session full
    figure; 
    subplot(2,2,1)
    hold on
    imagesc(all_common_cell_matrices{icomp,1}(sort_idx,:))
    set(gca,'TickLength',[0, 0]); box off;
    title([num2str(icomp) ' full'])
    ylabel('Neuron')
    xlabel('Time')
    
        % red line
        event_frame = cumsum([0 tses]).*100;
        xticks_hold = [];
        for i = event_frame 
           plot(i.*[1 1], ylim, 'r-', 'linewidth', 1.0) 
           xticks_hold = [xticks_hold i];
        end
        
        % aesthetics
        set(gca, 'YDir','reverse')
        ylim([0+0.5 size(all_common_cell_matrices{icomp,1}(sort_idx,:),1)+0.5])
        xlim([0.5 size(all_common_cell_matrices{icomp,1}(sort_idx,:),2)+0.5])
        xticks(xticks_hold)
        xticklabels_universe = {'NP1', 'NP2', 'HE1', 'TO', 'RD', 'RWD', 'TE'};
        xticklabels(xticklabels_universe)
    
    % plot first session zoom
    subplot(2,2,3)
    hold on
    imagesc(all_common_cell_matrices{icomp,1}(sort_idx, 1:zoom_portion(end)))
    set(gca,'TickLength',[0, 0]); box off;
    title([num2str(icomp) ' zoom'])
    ylabel('Neuron')
    xlabel('Time')
    
        % red line
        event_frame = cumsum([0 tses]).*100;
        xticks_hold = [];
        for i = event_frame 
           plot(i.*[1 1], ylim, 'r-', 'linewidth', 1.0) 
           xticks_hold = [xticks_hold i];
        end
        
        % aesthetics
        set(gca, 'YDir','reverse')
        ylim([0+0.5 size(all_common_cell_matrices{icomp,1}(sort_idx, 1:zoom_portion(end)),1)+0.5])
        xlim([0.5 size(all_common_cell_matrices{icomp,1}(sort_idx, 1:zoom_portion(end)),2)+0.5])
        xticks(xticks_hold)
        xticklabels_universe = {'NP1', 'NP2', 'HE1', 'TO', 'RD', 'RWD', 'TE'};
        xticklabels(xticklabels_universe)
        
    
    
    %sort and plot second session full
    subplot(2,2,2)
    hold on
    imagesc(all_common_cell_matrices{icomp,2}(sort_idx,:))
    set(gca,'TickLength',[0, 0]); box off;
    title([num2str(icomp+1) ' full'])
    ylabel('Neuron')
    xlabel('Time')
    
        % red line
        event_frame = cumsum([0 tses]).*100;
        xticks_hold = [];
        for i = event_frame 
           plot(i.*[1 1], ylim, 'r-', 'linewidth', 1.0) 
           xticks_hold = [xticks_hold i];
        end
        
        % aesthetics
        set(gca, 'YDir','reverse')
        ylim([0+0.5 size(all_common_cell_matrices{icomp,2}(sort_idx,:),1)+0.5])
        xlim([0.5 size(all_common_cell_matrices{icomp,2}(sort_idx,:),2)+0.5])
        xticks(xticks_hold)
        xticklabels_universe = {'NP1', 'NP2', 'HE1', 'TO', 'RD', 'RWD', 'TE'};
        xticklabels(xticklabels_universe)
    
    %sort and plot second session zoom
    subplot(2,2,4)
    hold on
    imagesc(all_common_cell_matrices{icomp,2}(sort_idx, 1:zoom_portion(end)))
    set(gca,'TickLength',[0, 0]); box off;
    title([num2str(icomp+1) ' zoom'])
    ylabel('Neuron')
    xlabel('Time')
    
        % red line
        event_frame = cumsum([0 tses]).*100;
        xticks_hold = [];
        for i = event_frame 
           plot(i.*[1 1], ylim, 'r-', 'linewidth', 1.0) 
           xticks_hold = [xticks_hold i];
        end
        
        % aesthetics
        set(gca, 'YDir','reverse')
        ylim([0+0.5 size(all_common_cell_matrices{icomp,2}(sort_idx, 1:zoom_portion(end)),1)+0.5])
        xlim([0.5 size(all_common_cell_matrices{icomp,2}(sort_idx, 1:zoom_portion(end)),2)+0.5])
        xticks(xticks_hold)
        xticklabels_universe = {'NP1', 'NP2', 'HE1', 'TO', 'RD', 'RWD', 'TE'};
        xticklabels(xticklabels_universe)

    set(gcf, 'Position', [219    18   894   863])
end






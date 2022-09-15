
%% details
mevar_subjs = [{'651049m1'} {'658648m2'} {'691359m1'} {'690330m1'} {'691359m2'} {'696944m3'} {'699437m4'} {'699438m3'}];
hivar_subjs = [{'687034m4'} {'683472m2'} {'699437m2'} {'699438m1'} {'683472m3'}];
num_mevar_subjs = length(mevar_subjs);
num_hivar_subjs = length(hivar_subjs);
probe_nums = [1:3:19 20]; 
green_color = [0 180 0]./255;
blue_color =  [0 115 207]./255;

%% compute or load cell correlations
%
%[all_cell_corrs_mevar, all_merge_mtx_mevar, subj_corr_cell_mevar, all_common_cell_matrices_mevar, subj_corr_mean{1s_mtx_mevar] = cell_turnover_timewarp_trials_multisubj('mevar', mevar_subjs);
%[all_cell_corrs_hivar, all_merge_mtx_hivar, subj_corr_cell_hivar, all_common_cell_matrices_hivar, subj_corr_means_mtx_hivar] = cell_turnover_timewarp_trials_multisubj('hivar', hivar_subjs);
%save('cell_corrs_all_may2022.mat', '-v7.3')

% OR
% load('cell_corrs_all_may2022.mat')
%}



%% cell level probe plots

% subject index
subj_idx_mevar = cell(size(all_cell_corrs_mevar));
for isubj = 1:num_mevar_subjs
    for isesh1 = 1:20
        for isesh2 = 1:20
            subj_idx_mevar{isesh1, isesh2} = [subj_idx_mevar{isesh1, isesh2}; repmat(isubj, size(subj_corr_cell_mevar{isubj}{isesh1, isesh2}))];
        end
    end
end
subj_idx_hivar = cell(size(all_cell_corrs_hivar));
for isubj = 1:num_hivar_subjs
    for isesh1 = 1:20
        for isesh2 = 1:20
            subj_idx_hivar{isesh1, isesh2} = [subj_idx_hivar{isesh1, isesh2}; repmat(isubj, size(subj_corr_cell_hivar{isubj}{isesh1, isesh2}))];
        end
    end
end

% get cell correlations for each probe session
probes_only_mevar = cell(8,8); 
probes_only_hivar = cell(8,8); 
for iprobe1 = 1:8
    for iprobe2 = 1:8
        probes_only_mevar{iprobe1, iprobe2} = all_cell_corrs_mevar{probe_nums(iprobe1), probe_nums(iprobe2)};
        probes_only_hivar{iprobe1, iprobe2} = all_cell_corrs_hivar{probe_nums(iprobe1), probe_nums(iprobe2)};
    end
end

% plot a mean correlation matrix
probes_only_mevar_mean_corm = cell2mat(cellfun(@mean,all_cell_corrs_mevar,'uni',0));
probes_only_hivar_mean_corm = cell2mat(cellfun(@mean,all_cell_corrs_hivar,'uni',0));
figure; imagesc(probes_only_mevar_mean_corm)
figure; imagesc(probes_only_hivar_mean_corm)

% get cell correlations for adjacent probe sessions
probes_only_adj_mevar = cell(1,8); probes_only_adj_mevar_subj_idx = cell(1,8);
probes_only_adj_hivar = cell(1,8); probes_only_adj_hivar_subj_idx = cell(1,8);
for iprobe_comp = 1:7
    probes_only_adj_mevar{iprobe_comp} = (probes_only_mevar{iprobe_comp, iprobe_comp+1}); 
    probes_only_adj_hivar{iprobe_comp} = (probes_only_hivar{iprobe_comp, iprobe_comp+1}); 
end

% plot adjacent probe cell corrs
figure; hold on
errorbar_plot_lineonly(probes_only_adj_mevar, 0, [], [], green_color)
errorbar_plot_lineonly(probes_only_adj_hivar, 0, [], [], blue_color)
for igroup_comp = 1:7
    [~,pval] = ttest2(probes_only_adj_mevar{igroup_comp}, probes_only_adj_hivar{igroup_comp});
    sig_asterisks(pval*7, igroup_comp-.1, .85)
end
tic_vect = [ -.9 -.8 -.6 -.3 0 .3 .6 .8 .9 ]; ylim(atanh([-.6 .9])); yticks(atanh(tic_vect)); yticklabels(tic_vect)
hold on; plot(xlim, [0 0], 'k--')
xlim([.5 7.5])
%}

% correlations with target probe
for target_probe = [1  8];
    probes_only_target_mevar = cell(num_mevar_subjs, 8);
    probes_only_target_hivar = cell(num_hivar_subjs, 8);
    for iprobe_comp = setdiff(1:8, target_probe)
        probes_only_target_mevar{iprobe_comp} = probes_only_mevar{target_probe, iprobe_comp};
        probes_only_target_hivar{iprobe_comp} = probes_only_hivar{target_probe, iprobe_comp};
    end

    % plot sessions from target
    figure; hold on;
    errorbar_plot_lineonly(probes_only_target_mevar, 0, [], [], green_color)
    errorbar_plot_lineonly(probes_only_target_hivar, 0, [], [], blue_color)
    plot(xlim, [0 0], 'k--')
    tic_vect = [ -.9 -.8 -.6 -.3 0 .3 .6 .8 .9 ]; ylim(atanh([-.6 .9])); yticks(atanh(tic_vect)); yticklabels(tic_vect)
    title(['Compare to probe ' num2str(target_probe)])
    for igroup_comp = setdiff(1:8, target_probe)
        [~,pval] = ttest2(probes_only_target_mevar{igroup_comp}, probes_only_target_hivar{igroup_comp})
        sig_asterisks(pval*7, igroup_comp-.1, .85)
    end
    xlabel('Probe')
    ylabel('Similarity')
end






%% subject level probe plots

% subject means must include at least min_cell_ct number of cells
min_cell_ct = 2;
for isesh1=1:20
    for isesh2=1:20
        for isubj_mevar = 1:num_mevar_subjs
            if length(subj_corr_cell_mevar{isubj_mevar}{isesh1, isesh2})<min_cell_ct
                subj_corr_cell_mevar{isubj_mevar}{isesh1, isesh2} = nan;
            end
            subj_corr_means_mtx_mevar(isesh1,isesh2,isubj_mevar) = mean(subj_corr_cell_mevar{isubj_mevar}{isesh1, isesh2});
        end
        for isubj_hivar = 1:num_hivar_subjs
            if length(subj_corr_cell_hivar{isubj_hivar}{isesh1, isesh2})<min_cell_ct
                subj_corr_cell_hivar{isubj_hivar}{isesh1, isesh2} = nan;
            end
            subj_corr_means_mtx_hivar(isesh1,isesh2,isubj_hivar) = mean(subj_corr_cell_hivar{isubj_hivar}{isesh1, isesh2});
        end
    end
end

% plot mean corms
figure; 
subplot(1,2,1); imagesc(nanmean(subj_corr_means_mtx_mevar,3)); axis square
set(gca,'TickLength',[0, 0]); xticks(probe_nums); yticks(probe_nums); 
caxis([-1 1]); title predictable; colorbar; xlabel('Session'); ylabel('Session')
subplot(1,2,2); imagesc(nanmean(subj_corr_means_mtx_hivar,3)); axis square
set(gca,'TickLength',[0, 0]); xticks(probe_nums); yticks(probe_nums);
caxis([-1 1]); title unpredictable; colorbar; xlabel('Session'); ylabel('Session')

% probes only
probes_only_subj_mevar = nan(8,8,num_mevar_subjs); 
probes_only_subj_hivar = nan(8,8,num_hivar_subjs); 
for iprobe1 = 1:8
    for iprobe2 = 1:8
        % extract
        probes_only_subj_mevar(iprobe1, iprobe2, min_cell_idx_mevar) = subj_corr_means_mtx_mevar(probe_nums(iprobe1), probe_nums(iprobe2), :);
        probes_only_subj_hivar(iprobe1, iprobe2, min_cell_idx_hivar) = subj_corr_means_mtx_hivar(probe_nums(iprobe1), probe_nums(iprobe2), :);
    end
end

% plot mean probes only corms
figure; 
subplot(1,2,1); imagesc(nanmean(probes_only_subj_mevar,3)); axis square
caxis([-1 1]);set(gca,'TickLength',[0, 0]); xlabel('Probe'); ylabel('Probe');
xticks(1:8); yticks(1:8); title predictable; colorbar
subplot(1,2,2); imagesc(nanmean(probes_only_subj_hivar,3)); axis square
caxis([-1 1]);set(gca,'TickLength',[0, 0]); xlabel('Probe'); ylabel('Probe');
xticks(1:8); yticks(1:8); title unpredictable; colorbar

% atahn
probes_only_subj_mevar = atanh(probes_only_subj_mevar);
probes_only_subj_hivar = atanh(probes_only_subj_hivar);

% adjacent
subj_probe_means_adjacent_mevar = nan(num_mevar_subjs, 7);
subj_probe_means_adjacent_hivar = nan(num_hivar_subjs, 7);
for iprobe_comp = 1:7
    subj_probe_means_adjacent_mevar(:,iprobe_comp) = probes_only_subj_mevar(iprobe_comp, iprobe_comp+1, :);
    subj_probe_means_adjacent_hivar(:,iprobe_comp) = probes_only_subj_hivar(iprobe_comp, iprobe_comp+1, :);
end

% plot adjacent
figure; hold on;
errorbar_mtx((subj_probe_means_adjacent_mevar), green_color)
errorbar_mtx((subj_probe_means_adjacent_hivar), blue_color)
plot(xlim, [0 0], 'k--')
tic_vect = [ -.9 -.8 -.6 -.3 0 .3 .6 .8 .9 ]; ylim(atanh([-.6 .9])); yticks(atanh(tic_vect)); yticklabels(tic_vect)
title('Adjacent')
for igroup_comp = 1:7
    [~,pval] = ttest2(subj_probe_means_adjacent_mevar(:, igroup_comp), subj_probe_means_adjacent_hivar(:, igroup_comp));
    sig_asterisks(pval, igroup_comp-.1, 1.15)
end
xlabel('Probe')
ylabel('Similarity')

% correlations with target probe
target_probe = 8;
subj_probe_means_target_mevar = nan(num_mevar_subjs, 8);
subj_probe_means_target_hivar = nan(num_hivar_subjs, 8);
for iprobe_comp = setdiff(1:8, target_probe)
    subj_probe_means_target_mevar(:,iprobe_comp) = probes_only_subj_mevar(target_probe, iprobe_comp, :);
    subj_probe_means_target_hivar(:,iprobe_comp) = probes_only_subj_hivar(target_probe, iprobe_comp, :);
end

% plot sessions from target
figure; hold on;
errorbar_mtx((subj_probe_means_target_mevar), green_color)
errorbar_mtx((subj_probe_means_target_hivar), blue_color)
plot(xlim, [0 0], 'k--')
tic_vect = [ -.9 -.8 -.6 -.3 0 .3 .6 .8 .9 ]; ylim(atanh([-.6 .9])); yticks(atanh(tic_vect)); yticklabels(tic_vect)
title(['Compare to probe ' num2str(target_probe)])
for igroup_comp = setdiff(1:8, target_probe)
    [~,pval] = ttest2(subj_probe_means_target_mevar(:, igroup_comp), subj_probe_means_target_hivar(:, igroup_comp));
    sig_asterisks(pval, igroup_comp-.1, 1.15)
end
xlabel('Probe')
ylabel('Similarity')

% correlations with target probe
target_probe = 1;
subj_probe_means_target_mevar = nan(num_mevar_subjs, 8);
subj_probe_means_target_hivar = nan(num_hivar_subjs, 8);
for iprobe_comp = setdiff(1:8, target_probe)
    subj_probe_means_target_mevar(:,iprobe_comp) = probes_only_subj_mevar(target_probe, iprobe_comp, :);
    subj_probe_means_target_hivar(:,iprobe_comp) = probes_only_subj_hivar(target_probe, iprobe_comp, :);
end

% plot sessions from target
figure; hold on;
errorbar_mtx((subj_probe_means_target_mevar), green_color)
errorbar_mtx((subj_probe_means_target_hivar), blue_color)
plot(xlim, [0 0], 'k--')
tic_vect = [ -.9 -.8 -.6 -.3 0 .3 .6 .8 .9 ]; ylim(atanh([-.6 .9])); yticks(atanh(tic_vect)); yticklabels(tic_vect)
title(['Compare to probe ' num2str(target_probe)])
for igroup_comp = setdiff(1:8, target_probe)
    [~,pval] = ttest2(subj_probe_means_target_mevar(:, igroup_comp), subj_probe_means_target_hivar(:, igroup_comp));
    sig_asterisks(pval, igroup_comp-.1, 1.15)
end
xlabel('Probe')
ylabel('Similarity')


%% subject level begining and end plots


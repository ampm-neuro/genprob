function probe_corr_mtx_twogroups(training_group1, training_group2)
% plut overlaid probe behavior for two groups, and a subplot of corms

num_freqs = 41;

ALL_probe_waits = cell(1,9);
ALL_subj_ids = cell(1,9);

for iprobe = 1:8
    figure; hold on
    [probe_waits_group1, subjects_ids_group1] = plot_probe_stages(iprobe, training_group1); 
    [probe_waits_group2, subjects_ids_group2] = plot_probe_stages(iprobe, training_group2); 
    
    title(['Probe 0' num2str(iprobe)]); 
    ylim([0 30]); 

    ALL_probe_waits_group1{iprobe} = reshape(probe_waits_group1, num_freqs, length(probe_waits_group1)/num_freqs)';
    ALL_subj_ids_group1{iprobe} = unique(subjects_ids_group1, 'rows' ,'stable');
    
    ALL_probe_waits_group2{iprobe} = reshape(probe_waits_group2, num_freqs, length(probe_waits_group2)/num_freqs)';
    ALL_subj_ids_group2{iprobe} = unique(subjects_ids_group2, 'rows' ,'stable');
end


% insert nans for missing subjects (uses first probe as master subj list)
for iprobe = 2:size(ALL_subj_ids_group1,2)
    
    %HERE
    nan_hold = nan(size(ALL_subj_ids_group1{1},1),num_freqs);
    nan_hold(ismember(ALL_subj_ids_group1{1}, ALL_subj_ids_group1{iprobe}, 'rows'), :) = ALL_probe_waits_group1{iprobe};
    ALL_probe_waits_group1{iprobe} = nan_hold;
    
end
for iprobe = 2:size(ALL_subj_ids_group2,2)
    
    nan_hold = nan(size(ALL_subj_ids_group2{1},1),num_freqs);
    nan_hold(ismember(ALL_subj_ids_group2{1}, ALL_subj_ids_group2{iprobe}, 'rows'), :) = ALL_probe_waits_group2{iprobe};
    ALL_probe_waits_group2{iprobe} = nan_hold;
    
end


% corr matrix for each subject
%{
corr_mtx = nan(size(ALL_probe_waits,2), size(ALL_probe_waits,2), size(ALL_subj_ids{1}, 1));

for isubj = 1:size(ALL_subj_ids{1},1)
    for iprobe1 = 1:size(ALL_probe_waits,2)        
        for iprobe2 = 1:size(ALL_probe_waits,2)

            waits_1 = ALL_probe_waits{iprobe1}(isubj,:)';
            waits_2 = ALL_probe_waits{iprobe2}(isubj,:)';
            
            nnan_idx = ~isnan(waits_1) & ~isnan(waits_2);
            
            if sum(nnan_idx)>10
                corr_mtx(iprobe1, iprobe2, isubj) = corr(waits_1(nnan_idx), waits_2(nnan_idx));
            end
            
        end
    end
end


figure; 
imagesc(nanmean(corr_mtx,3))
%}


% corr matrix for averages
corr_mtx_group1 = nan(size(ALL_probe_waits_group1,2), size(ALL_probe_waits_group1,2));
for iprobe1 = 1:size(ALL_probe_waits_group1,2)  
    for iprobe2 = 1:size(ALL_probe_waits_group1,2)
        waits_1 = nanmean(ALL_probe_waits_group1{iprobe1})';
        waits_2 = nanmean(ALL_probe_waits_group1{iprobe2})';
        nnan_idx = ~isnan(waits_1) & ~isnan(waits_2);
        if sum(nnan_idx)>10
            corr_mtx_group1(iprobe1, iprobe2) = corr(waits_1(nnan_idx), waits_2(nnan_idx));
        end
    end
end
corr_mtx_group2 = nan(size(ALL_probe_waits_group2,2), size(ALL_probe_waits_group2,2));
for iprobe1 = 1:size(ALL_probe_waits_group2,2)  
    for iprobe2 = 1:size(ALL_probe_waits_group2,2)
        waits_1 = nanmean(ALL_probe_waits_group2{iprobe1})';
        waits_2 = nanmean(ALL_probe_waits_group2{iprobe2})';
        nnan_idx = ~isnan(waits_1) & ~isnan(waits_2);
        if sum(nnan_idx)>10
            corr_mtx_group2(iprobe1, iprobe2) = corr(waits_1(nnan_idx), waits_2(nnan_idx));
        end
    end
end

% corm
figure; 
subplot(1,2,1)
ampm_pcolor(nanmean(corr_mtx_group1,3))
axis square
set(gca,'TickLength',[0, 0]);
caxis([-1 1])
ylabel('Probe number')
xlabel('Probe number')
title(training_group1)
colorbar

subplot(1,2,2)
ampm_pcolor(nanmean(corr_mtx_group2,3))
axis square
set(gca,'TickLength',[0, 0]);
caxis([-1 1])
ylabel('Probe number')
xlabel('Probe number')
title(training_group2)
colorbar



% corm fisher z transform
figure; 
subplot(1,2,1)
ampm_pcolor(atanh(nanmean(corr_mtx_group1,3)))
axis square
set(gca,'TickLength',[0, 0]);
caxis(atanh([-.925 .925]))
ylabel('Probe number')
xlabel('Probe number')
title(training_group1)
colorbar
subplot(1,2,2)
ampm_pcolor(atanh(nanmean(corr_mtx_group2,3)))
axis square
set(gca,'TickLength',[0, 0]);
caxis(atanh([-.925 .925]))
ylabel('Probe number')
xlabel('Probe number')
title(training_group2)
colorbar

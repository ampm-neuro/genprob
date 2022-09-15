function probe_corr_mtx_opto(training_group)
% separate correlation matrices of probes with laser ON and OFF


%training_group = 'train_mevar_optoExp_acc';

num_freqs = 41;

ALL_probe_waits_ON = cell(1,8);
ALL_probe_waits_OFF = cell(1,8);
ALL_subj_ids = cell(1,8);


% get all paths
for iprobe = 1:8
    if iprobe <=6
        path_str_1 = 'preprobe_opto';
        path_str_2 = ['_0' num2str(iprobe) '_'];
    elseif iprobe > 6
        path_str_1 = 'postprobe_opto';
        path_str_2 = ['_0' num2str(iprobe-6) '_'];
    end
    paths = get_file_paths_targeted('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\', training_group, path_str_1, path_str_2);
    
    
    % get wait times and subject IDs
    ALL_probe_waits_ON{iprobe} = nan(size(paths,1), num_freqs);
    ALL_probe_waits_OFF{iprobe} = nan(size(paths,1), num_freqs);
    for ipath = 1:size(paths,1)
    
        % load path
        load(paths{ipath}, 'trl_mtx')
        
        interp_method = 0; %no interp
        
        
        % subject id
        ALL_subj_ids{iprobe} = [ALL_subj_ids{iprobe}; find_subj_id(paths{ipath})];
        
        % laser ON
        ALL_probe_waits_ON{iprobe}(ipath,:) = wait_times_prep(trl_mtx(trl_mtx(:,13)==1,:), 2, 2)';

        % laser OFF
        ALL_probe_waits_OFF{iprobe}(ipath,:) = wait_times_prep(trl_mtx(trl_mtx(:,13)==0,:), 2, 2)';
        
        
        
    end
    
end

% insert nans for missing subjects (uses first probe as master subj list)
for iprobe = 2:size(ALL_subj_ids,2)
    
    % catch no subj
    if isempty(ALL_subj_ids) || isempty(ALL_subj_ids{iprobe})
        continue
    end
    
    % laser ON
    nan_hold = nan(size(ALL_subj_ids{1},1),num_freqs);
    nan_hold(ismember(ALL_subj_ids{1}, ALL_subj_ids{iprobe}, 'rows'), :) = ALL_probe_waits_ON{iprobe};
    ALL_probe_waits_ON{iprobe} = nan_hold;
    
    % laser OFF
    nan_hold = nan(size(ALL_subj_ids{1},1),num_freqs);
    nan_hold(ismember(ALL_subj_ids{1}, ALL_subj_ids{iprobe}, 'rows'), :) = ALL_probe_waits_OFF{iprobe};
    ALL_probe_waits_OFF{iprobe} = nan_hold;
end




% corr matrix for averages LASER ON
corr_mtx = nan(size(ALL_probe_waits_ON,2), size(ALL_probe_waits_ON,2));

for iprobe1 = 1:size(ALL_probe_waits_ON,2)  
    
    for iprobe2 = 1:size(ALL_probe_waits_ON,2)

        waits_1 = nanmean(ALL_probe_waits_ON{iprobe1})';
        waits_2 = nanmean(ALL_probe_waits_ON{iprobe2})';

        nnan_idx = ~isnan(waits_1) & ~isnan(waits_2);

        if sum(nnan_idx)>10
            corr_mtx(iprobe1, iprobe2) = corr(waits_1(nnan_idx), waits_2(nnan_idx));
        end

    end
end

figure; 
imagesc(nanmean(corr_mtx,3))
axis square
set(gca,'TickLength',[0, 0]);
caxis([-1 1])
ylabel('Probe number')
xlabel('Probe number')
title([training_group '; Laser ON'])


% corr matrix for averages LASER OFF
corr_mtx = nan(size(ALL_probe_waits_OFF,2), size(ALL_probe_waits_OFF,2));

for iprobe1 = 1:size(ALL_probe_waits_OFF,2)  
    
    for iprobe2 = 1:size(ALL_probe_waits_OFF,2)

        waits_1 = nanmean(ALL_probe_waits_OFF{iprobe1})';
        waits_2 = nanmean(ALL_probe_waits_OFF{iprobe2})';

        nnan_idx = ~isnan(waits_1) & ~isnan(waits_2);

        if sum(nnan_idx)>10
            corr_mtx(iprobe1, iprobe2) = corr(waits_1(nnan_idx), waits_2(nnan_idx));
        end

    end
end

figure; 
imagesc(nanmean(corr_mtx,3))
axis square
set(gca,'TickLength',[0, 0]);
caxis([-1 1])
ylabel('Probe number')
xlabel('Probe number')
title([training_group '; Laser OFF'])



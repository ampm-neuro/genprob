function ALL_audiodiscrim_zdiff
% plots zdiffs on last day of each audio discrim level

% path to audio discrimination behavior files
aud_discrim_path = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\audio_discrim';

% folders for each subject
subject_folders = dir(aud_discrim_path);
subject_folders(1:2) = [];

% number of folders
num_subjs = size(subject_folders,1);

% iteratively gather subject's files
subject_fpaths = cell(num_subjs,1);
for isubj = 1:num_subjs
    
    subj_folder_name = subject_folders(isubj).name;
    subject_fpaths{isubj} = get_file_paths_all...
        (['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\audio_discrim\' subj_folder_name]);

end

% identify discrimination levels
num_discrim_levels = 6;
subject_dlvl_fpaths = cell(num_subjs, num_discrim_levels);
for isubj = 1:num_subjs
    for idlvl = 1:num_discrim_levels
        subject_dlvl_fpaths{isubj,idlvl} = ...
            subject_fpaths{isubj}(contains(subject_fpaths{isubj}, ['DiscrimTest_0' num2str(idlvl-1)]));
    end
end

% final session at each discrimination level
subject_dlvl_fpaths_fin = cell(num_subjs, num_discrim_levels);
for isubj = 1:num_subjs
    for idlvl = 1:num_discrim_levels
        
        try
            subject_dlvl_fpaths_fin{isubj,idlvl} = subject_dlvl_fpaths{isubj,idlvl}{end};
        catch
        end
    end
end

% zdiffs
subj_dlvl_zdiffs = nan(num_subjs, num_discrim_levels);
for isubj = 1:num_subjs
    for idlvl = 1:num_discrim_levels
        
        % load session   
        if ~isempty(subject_dlvl_fpaths_fin{isubj,idlvl})
            load(subject_dlvl_fpaths_fin{isubj,idlvl}, 'trl_mtx')
        else
            continue
        end
        
        % find wait times
        [wait_times, frequencies] = wait_times_prep(trl_mtx, 2);
        rich_waits = wait_times(frequencies == min(frequencies));
        poor_waits = wait_times(frequencies ~= min(frequencies));
        
        % compute zifference
        subj_dlvl_zdiffs(isubj,idlvl) = zdiff(rich_waits,poor_waits);
        
    end
end


% plot
figure; hold on
errorbar_plot([{subj_dlvl_zdiffs(:,1)},{subj_dlvl_zdiffs(:,2)},{subj_dlvl_zdiffs(:,3)},{subj_dlvl_zdiffs(:,4)},{subj_dlvl_zdiffs(:,5)},{subj_dlvl_zdiffs(:,6)}]);

% aesthetics
set(gca,'TickLength',[0, 0]); box off;
xticks(1:6)
xticklabels({'0', '1', '2', '3', '4', '5'})
xlabel('Units between rich and poor tones')
ylim([-1 5])
yticks(0:5)
hold on; plot(xlim, [1 1].*0, 'k--')
ylabel('Rich - Poor wait times (z)')




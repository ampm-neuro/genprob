function [coefEsts_out_normal, coefEsts_out_log] = probe_coefs(folderpath)
% iterate through subjects and training stages to determine how similar
% performance on the first test problem of a training stage is to
% performance on the probe session the day before



%% get unique subjects
folder_contents = dir(folderpath);
folder_contents = struct2cell(folder_contents);
unq_subjects = folder_contents(1,contains(folder_contents(1,:), 'm'));

%% get session files

% preallocate cells to hold session file paths
paths_probes = cell(length(unq_subjects),1);

% iterate through subjects
for isubj = 1:length(unq_subjects)
    paths_probes{isubj} = get_file_paths_targeted(folderpath, {'preprobe', unq_subjects{isubj}});
    paths_probes{isubj} = [paths_probes{isubj}; get_file_paths_targeted(folderpath, {'postprobe_01', unq_subjects{isubj}})];
end



%% Compute ceofs for each probe
coefEsts_out_normal = nan(7, length(unq_subjects));
coefEsts_out_log = nan(7, length(unq_subjects));
for isubj = 1:length(paths_probes)
    for iprobe = 1:size(paths_probes{isubj},1)
    
        % load probe session
        load(paths_probes{isubj}{iprobe})

        % compute coefficients
        [~, ~, ~, ~, coefEsts] = wait_times_prep(trl_mtx, 2, 0);

        % load coefficient outputs
        coefEsts_out_normal(iprobe,isubj) = coefEsts(4); % normal
        coefEsts_out_log(iprobe,isubj) = coefEsts(3); % log
    
    end
end
            





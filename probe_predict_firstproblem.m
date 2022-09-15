function [train_mean_waits, probe_model_waits] = probe_predict_firstproblem(folderpath)
% iterate through subjects and training stages to determine how similar
% performance on the first test problem of a training stage is to
% performance on the probe session the day before


%% hardcoded parameters

%hard code frequencies 
hcf = [5000,5249,5511,5786,6074,6377,6695,7028,7379,7747,8133,8538,...
    8964,9411,9880,10372,10890,11432,12002,12601,13229,13888,14581,...
    15307,16070,16872,17713,18596,19523,20496,21518,22590,23716,24899,...
    26140,27443,28811,30247,31755,33338,35000];


%% get unique subjects
folder_contents = dir(folderpath);
folder_contents = struct2cell(folder_contents);
unq_subjects = folder_contents(1,contains(folder_contents(1,:), 'm'));

%% get session files

% preallocate cells to hold session file paths
paths_preprobes_I = cell(length(unq_subjects),1);
paths_train01s_I = cell(length(unq_subjects),1);

% iterate through subjects
for isubj = 1:length(unq_subjects)
    % pre probe session files
    %get_file_paths_targeted(folderpath, {'preprobe', unq_subjects{isubj}})
    paths_preprobes_I{isubj} = get_file_paths_targeted(folderpath, {'preprobe', unq_subjects{isubj}});
    % first training session files
    %get_file_paths_targeted(folderpath, {'_train', '-01', unq_subjects{isubj}})
    train_paths_mevar = get_file_paths_targeted(folderpath, {'_mevar', '-01', unq_subjects{isubj}});
    train_paths_hivar = get_file_paths_targeted(folderpath, {'_hivar', '-01', unq_subjects{isubj}});
    if size(train_paths_mevar,1) > size(train_paths_hivar,1)
        paths_train01s_I{isubj} = train_paths_mevar;
    else
        paths_train01s_I{isubj} = train_paths_hivar;    
    end
    paths_train01s_I{isubj} = paths_train01s_I{isubj}(~contains(paths_train01s_I{isubj}, 'notone'));
end


%% filter session files to only include matching pre-probes and training sessions

% iterate though subjects
max_train_sesh = nan(length(unq_subjects),1);
for isubj = 1:length(unq_subjects)
    
    % find maximum training session
    train_sesh_nums = [];
    for itrain = 1:length(paths_train01s_I{isubj})
                
        if size(train_paths_mevar,1) > size(train_paths_hivar,1)
            str_find_str = '_mevar';
        else
            str_find_str = '_hivar';
        end
            
        str_idx = strfind(paths_train01s_I{isubj}{itrain}, str_find_str);
        str_idx = str_idx(end);
        str_idx = str_idx + length(str_find_str);
        train_sesh_num = paths_train01s_I{isubj}{itrain}(str_idx : str_idx+1);
        train_sesh_num = str2num(train_sesh_num);
        train_sesh_nums = [train_sesh_nums; train_sesh_num];
    end    
    max_train_sesh(isubj) = max(train_sesh_nums);
    
    % preallocate cell to hold session file paths
    paths_preprobes_II = cell(max_train_sesh(isubj), 1);
    paths_train01s_II = cell(max_train_sesh(isubj), 1);
    
    % iterate through potential training sessions
    for itrain = 1:max_train_sesh(isubj)
        
        % training session number string
        train_str = num2str(itrain);
        if length(train_str)<2
            train_str = ['0' train_str];
        end

        % skip subject if there are no probe or training sessions
        if isempty(paths_preprobes_I{isubj}) || isempty(paths_train01s_I{isubj})
            continue
        end
        
        % load only paired session files
        probe_file = paths_preprobes_I{isubj}(contains(paths_preprobes_I{isubj}, ['_preprobe_' train_str]));
        if size(train_paths_mevar,1) > size(train_paths_hivar,1)
            train_file = paths_train01s_I{isubj}(contains(paths_train01s_I{isubj}, ['_mevar' train_str]));
        else
            train_file = paths_train01s_I{isubj}(contains(paths_train01s_I{isubj}, ['_hivar' train_str]));
        end
        if length(probe_file)==1 && length(train_file) == 1
            paths_preprobes_II{itrain} = probe_file{:};
            paths_train01s_II{itrain} = train_file{:};
        end
    end
    
    % load back into original cells
    paths_preprobes_I{isubj} = paths_preprobes_II;
    paths_train01s_I{isubj} = paths_train01s_II;

end


%% Compute correlation between probe and training wait times

% preallocate
train_mean_waits = cell(length(unq_subjects),1);
probe_model_waits = cell(length(unq_subjects),1);

% rich and poor hi-low order 1==rich, 2==poor
rich_poor_order = rich_bounds_prob('mevar', 0);
for irow = 1:size(rich_poor_order,1)
    if rich_poor_order(irow,1)<rich_poor_order(irow,2)
        rich_poor_order(irow,:) = [1 2];
    else
        rich_poor_order(irow,:) = [2 1];
    end
end

% iterate though subjects
for isubj = 1:length(unq_subjects)
    probe_model_waits_hold = [];
    train_mean_waits_hold = [];
    
    % iterate through session pairs
    for isesh = 1:size(paths_train01s_I{isubj},1)
    
        % load training session or skip
        %
        if ~isempty(paths_train01s_I{isubj}{isesh})
            load(paths_train01s_I{isubj}{isesh})
        else
            continue
        end            
            
        
        % compute mean training wait times
        %
        load(paths_train01s_I{isubj}{isesh})
        [mean_train_wait_times, ~, mean_train_tone_nums] = wait_times_prep(trl_mtx, 1, 4);
        
        
        %compute mean probe wait times
        %
        load(paths_preprobes_I{isubj}{isesh})
        [~, ~, ~, modelFun, coefEsts] = wait_times_prep(trl_mtx, 2, 0);

            % model wait times at frequencies seen in training session   
            modeled_wait_times = modelFun(coefEsts, unique(mean_train_tone_nums));
        
            
            % MAKE SURE FREQUENCIES MATCH UP
            % yes, because modeled_wait_times were computed with training_
            % wait_tone nums as input
            
        % load
        %
        modeled_wait_times = modeled_wait_times(~isnan(mean_train_wait_times))';    
        mean_train_wait_times = mean_train_wait_times(~isnan(mean_train_wait_times))';

        % arrange by rich and poor
        modeled_wait_times = modeled_wait_times(rich_poor_order(isesh,:));
        mean_train_wait_times = mean_train_wait_times(rich_poor_order(isesh,:));
        
        % load output
        probe_model_waits_hold = [probe_model_waits_hold; modeled_wait_times];
        train_mean_waits_hold = [train_mean_waits_hold; mean_train_wait_times];
        
    end
    
    % load subject
    probe_model_waits{isubj} = probe_model_waits_hold;
    train_mean_waits{isubj} = train_mean_waits_hold;
    
end


%% plot shit
%{
% preallocate
training_stage_probe_waits = cell(nanmax(max_train_sesh),1);
training_stage_train_waits = cell(nanmax(max_train_sesh),1);

% iterate through subjects
for isubj = 1:length(unq_subjects) 
    % iterate through training stages
    for itrain = 1:length(train_mean_waits{isubj})        
        % gather probe and training wait times
        if ~isempty(probe_model_waits{isubj}{itrain}) && ~isempty(train_mean_waits{isubj}{itrain})
            training_stage_probe_waits{itrain} = [training_stage_probe_waits{itrain}; probe_model_waits{isubj}{itrain}(:,2)]
            training_stage_train_waits{itrain} = [training_stage_train_waits{itrain}; train_mean_waits{isubj}{itrain}(:,2)]
        end
    end
end

% plot fit lines for each training stage
figure; hold on
for itrain = 1:length(training_stage_probe_waits)
    subplot(1,length(training_stage_probe_waits), itrain)

    try
    [R, p, se] = fit_line(training_stage_probe_waits{itrain}, training_stage_train_waits{itrain});
    axis([-2 5 -2 5])
    axis square

    title(['r=' num2str(R) '; p=' num2str(p) '; mse=' num2str(mean(se))])
    catch
    end
    
end
%}





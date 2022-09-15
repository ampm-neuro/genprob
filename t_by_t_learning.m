%[suprise_cell, discrim_cell, suprise_mtx_interp, problem_paths] = t_by_t_learning

% plot rich-tone preference overlaid with suprise
subject_ids = [unique_subjs('\train_mevar\');unique_subjs('\train_hivar\')] ;
%subject_ids = unique_subjs('\train_mevar\');
%subject_ids = unique_subjs('\train_mevar_imaging_hpc\');

%% basic info
wdw_size = [4 4];
load('unqfrq41', 'unqfrq41')

% behavioral folder
beh_fld = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\';

% get all first-day paths for each subject
problem_paths = cell(1, size(subject_ids,1));
for isubj = 1:size(subject_ids,1)
    all_subj_paths = get_file_paths_targeted(beh_fld, {subject_ids(isubj,:)});
    for iproblem = 1:6
        all_problem_files = [all_subj_paths(contains(all_subj_paths, ['var0' num2str(iproblem)])); ...
            all_subj_paths(contains(all_subj_paths, ['ctl0' num2str(iproblem)]))];
        problem_paths{isubj} = [problem_paths{isubj}; all_problem_files(1)];
    end
end
clearvars all_problem_files all_subj_paths iproblem isubj
  
    
%% iterate through learning sessions getting instantaneous rich-tone preferences
all_subj_wait_prefs = cell(1, size(subject_ids,1));
all_subj_problem_idx = cell(1, size(subject_ids,1));
for isubj = 1:length(problem_paths)
    
    % iterate through first-day sessions
    for isesh = 1:size(problem_paths{isubj},1)
        % compute trial-by-trial preference for rich tone
        [wait_prefs_sesh] = continuous_learning_singlesesh_suprise(problem_paths{isubj}{isesh}, wdw_size);
        all_subj_wait_prefs{isubj}{isesh} = wait_prefs_sesh;
        all_subj_problem_idx{isubj}{isesh} = repmat(isesh, size(wait_prefs_sesh));
    end
    
end


%% plot instantaneous preferences
all_subj_wait_prefs_interp = all_subj_wait_prefs;
all_pref_mtx = [];
for isubj = 1:length(problem_paths)
    for isesh = 1:length(problem_paths{isubj})
        % warp preferences
        preference_vect = all_subj_wait_prefs_interp{isubj}{isesh};
        all_subj_wait_prefs_interp{isubj}{isesh} = interp1(1:length(preference_vect), preference_vect, linspace(1,length(preference_vect),100));
    end
    
    % cell2mat
    all_pref_mtx = [all_pref_mtx; cell2mat(all_subj_wait_prefs_interp{isubj})];
end

% means and stds
all_pref_mean = nanmean(all_pref_mtx);
all_pref_ste = nanstd(all_pref_mtx)./sqrt(sum(~isnan(all_pref_mtx)));

%
figure; hold on
for isubj = 1:size(all_pref_mtx,1)
    plot(all_pref_mtx(isubj,:), 'color', .9.*[1 1 1], 'linewidth', 1)
end
plot(all_pref_mean, 'color', .4.*[1 1 1], 'linewidth', 2)
    plot(all_pref_mean-all_pref_ste, 'color', .65.*[1 1 1], 'linewidth', 1)
    plot(all_pref_mean+all_pref_ste, 'color', .65.*[1 1 1], 'linewidth', 1)  

plot(xlim, [0 0], 'k--');
ylim([-1 1]); set(gca,'TickLength',[0, 0]);
plot([100.5:100:600]'.*[1 1], ylim, 'k--')
title('Instantaneous preference for rich tone')


%% get pre-probe paths to identify prediction

% get all probe paths for each subject
probe_paths = cell(1, size(subject_ids,1));
for isubj = 1:size(subject_ids,1)
    probe_paths{isubj} = get_file_paths_targeted(beh_fld, {subject_ids(isubj,:), 'preprobe'});
end


%% iterate through probe sessions modelling future rich-tone preferences

% rich tones
rich_tones_predict = rich_bounds_prob('mevar', 0, 1);
rich_tones_unpredict = rich_bounds_prob('hivar', 0, 1);

% iterate through subjects
all_subj_wait_predictions = cell(1, size(subject_ids,1));
all_subj_preprobe_idx = cell(1, size(subject_ids,1));
for isubj = 1:length(probe_paths)
    % iterate through probe sessions
    for isesh = 1:size(probe_paths{isubj},1)
        % load session
        load(probe_paths{isubj}{isesh}, 'trl_mtx')
        % observered preprobe wait times
        [wait_times, ~, obs_freq_nums] = wait_times_prep(trl_mtx, 1, 0);
            try
                [~, coefEsts, modelFun] = ampm_normal_logistic_fit_algo(obs_freq_nums, wait_times, [nanmean(wait_times) 21 0 1]);
            catch
                try
                    [~, coefEsts, modelFun] = ampm_normal_logistic_fit_algo(obs_freq_nums, wait_times, [nanmean(wait_times) 20.5 0 1]);
                catch
                    try
                        [~, coefEsts, modelFun] = ampm_normal_logistic_fit_algo(obs_freq_nums, wait_times, [nanmean(wait_times) 21.5 0 1]);
                    catch
                        try
                            [~, coefEsts, modelFun] = ampm_normal_logistic_fit_algo(obs_freq_nums, wait_times, [nanmean(wait_times) 20.0 0 1]);
                        catch
                            try
                                [~, coefEsts, modelFun] = ampm_normal_logistic_fit_algo(obs_freq_nums, wait_times, [nanmean(wait_times) 22.0 0 1]);
                            catch
                            end
                        end
                    end
                end
            end
        % rich and poor predictions
        if contains(probe_paths{isubj}{isesh}, 'mevar')
            rich_poor_tones = rich_tones_predict(isesh,:);
        else
            rich_poor_tones = rich_tones_unpredict(isesh,:);
        end
        rich_poor_predictions = modelFun(coefEsts, rich_poor_tones);
        all_subj_wait_predictions{isubj}{isesh} = (rich_poor_predictions(1)-rich_poor_predictions(2))/(rich_poor_predictions(1)+rich_poor_predictions(2));
        all_subj_problem_idx{isubj}{isesh} = isesh;
    end
end



%% plot similarity to prediction
all_subj_sim_to_prediction = all_subj_wait_prefs;
all_subj_sim_to_prediction_interp = all_subj_wait_prefs_interp;
all_subj_sim_to_prediction_delta = all_subj_wait_prefs;
all_subj_sim_to_prediction_delta_interp = all_subj_wait_prefs_interp;
all_sim2pred_mtx = [];
all_sim2pred_mtx_delta = [];
for isubj = 1:length(problem_paths)
    for isesh = 1:length(problem_paths{isubj})
        % compute instantaneous difference from prediction
        all_subj_sim_to_prediction{isubj}{isesh} = all_subj_sim_to_prediction{isubj}{isesh} - all_subj_wait_predictions{isubj}{isesh};
        all_subj_sim_to_prediction_interp{isubj}{isesh} = all_subj_sim_to_prediction_interp{isubj}{isesh} - all_subj_wait_predictions{isubj}{isesh};
        all_subj_sim_to_prediction_delta{isubj}{isesh} = all_subj_sim_to_prediction{isubj}{isesh}(2:end) - all_subj_sim_to_prediction{isubj}{isesh}(1:end-1);
        all_subj_sim_to_prediction_delta_interp{isubj}{isesh} = all_subj_sim_to_prediction_interp{isubj}{isesh}(2:end) - all_subj_sim_to_prediction_interp{isubj}{isesh}(1:end-1);
    end
    
    % cell2mat
    all_sim2pred_mtx = [all_sim2pred_mtx; cell2mat(all_subj_sim_to_prediction_interp{isubj})];
    all_sim2pred_mtx_delta = [all_sim2pred_mtx_delta; cell2mat(all_subj_sim_to_prediction_delta_interp{isubj})];
end

% means and stds
all_sim2pred_mean = nanmean(all_sim2pred_mtx);
all_sim2pred_ste = nanstd(all_sim2pred_mtx)./sqrt(sum(~isnan(all_sim2pred_mtx)));

%
figure; hold on
for isubj = 1:size(all_sim2pred_mtx,1)
    plot(all_sim2pred_mtx(isubj,:), 'color', .9.*[1 1 1], 'linewidth', 1)
end
plot(all_sim2pred_mean, 'color', .4.*[1 1 1], 'linewidth', 2)
    plot(all_sim2pred_mean-all_sim2pred_ste, 'color', .65.*[1 1 1], 'linewidth', 1)
    plot(all_sim2pred_mean+all_sim2pred_ste, 'color', .65.*[1 1 1], 'linewidth', 1)  

plot(xlim, [0 0], 'k--');
ylim([-1 1]); set(gca,'TickLength',[0, 0]);
plot([100.5:100:600]'.*[1 1], ylim, 'k--')
title('Instantaneous preference for rich tone - prediction')

% means and stds
%
all_sim2pred_delta_mean = nanmean(all_sim2pred_mtx_delta);
all_sim2pred_delta_ste = nanstd(all_sim2pred_mtx_delta)./sqrt(sum(~isnan(all_sim2pred_mtx_delta)));

figure; hold on
for isubj = 1:size(all_sim2pred_mtx_delta,1)
    plot(all_sim2pred_mtx_delta(isubj,:), 'color', .9.*[1 1 1], 'linewidth', 1)
end
plot(all_sim2pred_delta_mean, 'color', .4.*[1 1 1], 'linewidth', 2)
    plot(all_sim2pred_delta_mean-all_sim2pred_delta_ste, 'color', .65.*[1 1 1], 'linewidth', 1)
    plot(all_sim2pred_delta_mean+all_sim2pred_delta_ste, 'color', .65.*[1 1 1], 'linewidth', 1)  

plot(xlim, [0 0], 'k--');
ylim([-1 1]); set(gca,'TickLength',[0, 0]);
plot([100.5:100:600]'.*[1 1], ylim, 'k--')
title('Instantaneous preference for rich tone - prediction (DELTA)')
%}



%% iterate through learning sessions identifying suprise
all_subj_suprise = cell(1, size(subject_ids,1));
all_subj_suprise_inst = cell(1, size(subject_ids,1));
for isubj = 1:length(problem_paths)
    
    % iterate through first-day sessions
    for isesh = 1:size(problem_paths{isubj},1)
        
        % iterate through trials
        load(problem_paths{isubj}{isesh}, 'trl_mtx')
        all_subj_suprise{isubj}{isesh} = nan(size(trl_mtx,1),1);
        all_subj_suprise_inst{isubj}{isesh} = nan(size(trl_mtx,1),1);
        
        % rich tone
        rich_tone = mode(floor(trl_mtx(trl_mtx(:,3)==1,2)));
        rich_idx = floor(trl_mtx(:,2))==rich_tone;
        %reward idx
        rwd_idx = ~isnan(trl_mtx(:,11));
        
        % suprise idx
        all_subj_suprise{isubj}{isesh}(rich_idx & rwd_idx) = 1; % yes rich, yes rwd
        all_subj_suprise{isubj}{isesh}(~rich_idx & ~rwd_idx) = 2; % no rich, no rwd
        all_subj_suprise{isubj}{isesh}(rich_idx & ~rwd_idx) = 3; % yes rich, no rwd
        all_subj_suprise{isubj}{isesh}(~rich_idx & rwd_idx) = 4; % no rich, yes rwd
        
        % compute proportion of suprising trials in window
        %{
        for itrl = 1 : size(trl_mtx,1)
            local_window = (itrl-wdw_size(1) : itrl+wdw_size(2)); % all past trials 
            %local_window = (1 : itrl); % all past trials 
            %local_window = (itrl+1:size(trl_mtx,1)); % future trials (not sig)
            local_window = local_window(local_window>0 & local_window<=size(trl_mtx,1));
            
            unpredicted_events = sum(ismember(all_subj_suprise{isubj}{isesh}(local_window), [3 4]));
            predicted_events = sum(ismember(all_subj_suprise{isubj}{isesh}(local_window), [1 2]));
            
            all_subj_suprise_inst{isubj}{isesh}(itrl) = (unpredicted_events-predicted_events)/(unpredicted_events+predicted_events);
        end
        
        % multiply suprise by subjects instantaneous prediction
        %all_subj_suprise_inst_fin{isubj}{isesh} = all_subj_suprise_inst{isubj}{isesh};
        all_subj_suprise_inst_fin{isubj}{isesh} = all_subj_suprise_inst{isubj}{isesh}.*all_subj_wait_prefs{isubj}{isesh};
                %all_subj_suprise_inst_fin{isubj}{isesh} = all_subj_suprise_inst{isubj}{isesh}.*(all_subj_wait_prefs{isubj}{isesh}./abs(all_subj_wait_prefs{isubj}{isesh}));
        %}
    end
    
end


%% Iterate through sessions identifying trial wait times
all_subj_waits = cell(1, size(subject_ids,1));
for isubj = 1:length(problem_paths)
    
    % iterate through first-day sessions
    for isesh = 1:size(problem_paths{isubj},1)
        
        % iterate through trials
        load(problem_paths{isubj}{isesh}, 'trl_mtx')
        all_subj_waits{isubj}{isesh} = trl_mtx(:,12);

    end
end



%% GLM of wait times and suprise trial type

instpref_rel2pred_cell = cell(1,6);
instpref_rel2pred_delta_cell = cell(1,6);
wait_times_cell = cell(1,6); 
suprise_cell = cell(1,6); 
subject_cell = cell(1,6); 
session_cell = cell(1,6);

for iprobe = 1:length(all_subj_suprise{1})
    % predict
    for isubj = 1:length(all_subj_suprise)
        % data
        instpref_rel2pred_cell{1,iprobe} = [instpref_rel2pred_cell{1,iprobe}; all_subj_sim_to_prediction{isubj}{iprobe}];
        instpref_rel2pred_delta_cell{1,iprobe} = [instpref_rel2pred_delta_cell{1,iprobe}; all_subj_sim_to_prediction_delta{isubj}{iprobe}];
        wait_times_cell{1,iprobe} = [wait_times_cell{1,iprobe}; all_subj_waits{isubj}{iprobe}(1:end-1)];
        
        suprise_windowct = running_window(all_subj_suprise{isubj}{iprobe}, [8,0]) ;
        suprise_cell{1,iprobe} = [suprise_cell{1,iprobe}; suprise_windowct(1:end-1)];
        % indices
        local_size = size(all_subj_suprise{isubj}{iprobe}(1:end-1));
        subject_cell{1,iprobe} = [subject_cell{1,iprobe}; repmat(isubj, local_size)];
        session_cell{1,iprobe} = [session_cell{1,iprobe}; repmat(iprobe, local_size)];
    end
end

% cell2mat
instpref_rel2pred = cell2mat(instpref_rel2pred_cell');
instpref_rel2pred_delta = cell2mat(instpref_rel2pred_delta_cell');

wait_times = cell2mat(wait_times_cell');
suprise = cell2mat(suprise_cell');
subject = cell2mat(subject_cell');
session = cell2mat(session_cell');




model_str = 'RwdPrefDelta~SupriseType + (1|Subject) + (1|Session)';
tbl = table(instpref_rel2pred_delta, wait_times, suprise, subject, session, 'VariableNames', {'RwdPrefDelta', 'WaitTimes', 'SupriseType', 'Subject', 'Session'});
tbl.Subject = categorical(tbl.Subject);
%tbl.SupriseType = categorical(tbl.SupriseType);
%fit = fitlm(tbl,model_str)
fit = fitlme(tbl, model_str)



%% internal function
function outvect = running_window(invect, windowsize)
     outvect = nan(size(invect));
     for it = 1:length(invect)
         wdw = (it-windowsize(1)):(it+windowsize(2));
         wdw = wdw(wdw>=1 & wdw<=length(invect));
         if ~isempty(wdw)
            outvect(it) = (sum(ismember(invect(wdw), [3 4])) - sum(ismember(invect(wdw), [1 2]))) / (sum(ismember(invect(wdw), [3 4])) + sum(ismember(invect(wdw), [1 2])));
         else
            outvect(it) = 0;
         end
     end

end











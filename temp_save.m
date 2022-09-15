

%% correct for time cell level
%{
% time between imaging session (in days)
%{
days_mtx_predict = [];
for isubj_predict = 1:length(predict_subjs)
    days_mtx_predict = cat(3, days_mtx_predict, session_days(predict_subjs{isubj_predict}));
end
days_mtx_unpredict = [];
for isubj_unpredict = 1:length(unpredict_subjs)
    days_mtx_unpredict = cat(3, days_mtx_unpredict, session_days(unpredict_subjs{isubj_unpredict}));
end
figure; subplot(1,2,1); imagesc(nanmean(days_mtx_predict,3)); caxis([0 60]);...
 subplot(1,2,2); imagesc(nanmean(days_mtx_unpredict,3)); caxis([0 60]);
save('days_between.mat', 'days_mtx_predict', 'days_mtx_unpredict')
%}
%
load('days_between.mat', 'days_mtx_predict', 'days_mtx_unpredict')
load('firing_field_data_predict')
load('firing_field_data_unpredict')
load('sub_corr_cells_atanh', 'subj_corr_cell_predict', 'subj_corr_cell_unpredict')

% remove identity diaganol
subj_corr_means_mtx_predict(subj_corr_means_mtx_predict>0.99) = nan; 
subj_corr_means_mtx_unpredict(subj_corr_means_mtx_unpredict>0.99) = nan; 
days_mtx_predict = days_mtx_predict(1:19,1:19,:);
    days_mtx_predict(days_mtx_predict==0) = nan;
days_mtx_unpredict = days_mtx_unpredict(1:19,1:19,:);
    days_mtx_unpredict(days_mtx_unpredict==0) = nan;
    
% assign a days-between value to each r value
predict_rvals_and_days = [];
for isubj = 1:length(subj_corr_cell_predict)
    for isesh1 = 1:size(subj_corr_cell_predict{isubj},1)
        for isesh2 = 1:size(subj_corr_cell_predict{isubj},2)
            % skip redundancies
            if isesh2<=isesh1
                continue
            end
            predict_rvals_and_days = [predict_rvals_and_days; subj_corr_cell_predict{isubj}{isesh1, isesh2} repmat(days_mtx_predict(isesh1, isesh2, isubj), size(subj_corr_cell_predict{isubj}{isesh1, isesh2}))];
        end
    end
end
unpredict_rvals_and_days = [];
for isubj = 1:length(subj_corr_cell_unpredict)
    for isesh1 = 1:size(subj_corr_cell_unpredict{isubj},1)
        for isesh2 = 1:size(subj_corr_cell_unpredict{isubj},2)
            % skip redundancies
            if isesh2<=isesh1
                continue
            end
            unpredict_rvals_and_days = [unpredict_rvals_and_days; subj_corr_cell_unpredict{isubj}{isesh1, isesh2} repmat(days_mtx_unpredict(isesh1, isesh2, isubj), size(subj_corr_cell_unpredict{isubj}{isesh1, isesh2}))];
        end
    end
end    
    

% combine
all_rvals_and_days = [predict_rvals_and_days; unpredict_rvals_and_days];
%all_rvals_and_days = unpredict_rvals_and_days;
    
% model time
ibins = [0:0.25:5 5.5:.5:30 31:1:40 42:2:60 65:5:100];
ibin_cell = cell(1,length(ibins)-1);
for ibin = 1:length(ibin_cell)
    ibin_cell{ibin} = all_rvals_and_days(all_rvals_and_days(:,2)>ibins(ibin) & all_rvals_and_days(:,2)<=ibins(ibin+1), 1);
end
max_x = 100;
[~, coefEsts, modelFun] = ampm_normal_fit(all_rvals_and_days(~isnan(all_rvals_and_days(:,2)) & all_rvals_and_days(:,2)<max_x, 2), all_rvals_and_days(~isnan(all_rvals_and_days(:,2)) & all_rvals_and_days(:,2)<max_x, 1), [0 1 0 150]);
coefEsts
% correct population overlap matrices
subj_corr_means_mtx_predict = subj_corr_means_mtx_predict - modelFun(coefEsts, days_mtx_predict);
subj_corr_means_mtx_unpredict = subj_corr_means_mtx_unpredict - modelFun(coefEsts, days_mtx_unpredict);
%

%plot data points, boxcar mean, and fit curve
%
figure; hold on; 
plot(all_rvals_and_days(:,2), all_rvals_and_days(:,1), 'o', 'color', 0.8.*[1 1 1])
errorbar_plot_lineonly(ibin_cell, 0, ibins(1:end-1), 0.8.*[1 1 1])
plot(ibins, modelFun(coefEsts, ibins), 'r-')
xlim([-1 90]); plot(xlim, [0 0], 'k--');
ylabel('Activity correlation')
xlabel('Days between imaging sessions')
tic_vect = [-.999 -.99 -.9 -.8 -.6 -.3 0 .3 .6 .8 .9 .99 .999]; ylim(atanh([-.999 .999])); yticks(atanh(tic_vect)); yticklabels(tic_vect)
%var_name = 'Time_curve_corr'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
%

% correct subject cells
%
for isubj = 1:length(subj_corr_cell_predict)
    for isesh1 = 1:size(subj_corr_cell_predict{isubj},1)
        for isesh2 = 1:size(subj_corr_cell_predict{isubj},2)
            subj_corr_cell_predict{isubj}{isesh1, isesh2} = subj_corr_cell_predict{isubj}{isesh1, isesh2} - modelFun(coefEsts, days_mtx_predict(isesh1, isesh2, isubj));
        end
    end
end
for isubj = 1:length(subj_corr_cell_unpredict)
    for isesh1 = 1:size(subj_corr_cell_unpredict{isubj},1)
        for isesh2 = 1:size(subj_corr_cell_unpredict{isubj},2)
            subj_corr_cell_unpredict{isubj}{isesh1, isesh2} = subj_corr_cell_unpredict{isubj}{isesh1, isesh2} - modelFun(coefEsts, days_mtx_unpredict(isesh1, isesh2, isubj));
        end
    end
end
save('sub_corr_cells_timecorrect', 'subj_corr_cell_predict', 'subj_corr_cell_unpredict')
%}


%% Cell level stage errorbar plots
%{

% load all held cells corrs
load('firing_field_data_predict')
load('firing_field_data_unpredict')
load('sub_corr_cells_timecorrect', 'subj_corr_cell_predict', 'subj_corr_cell_unpredict')

% overall neuron ids 
base_num = 0;
for isubj_predict = 1:length(predict_subjs)
    
    % unique cell numbers
    unique_numbers = unique(cell2mat(subj_corr_cell_cellid_predict{isubj_predict}(:)));
    unique_numbers = [unique_numbers (1:length(unique_numbers))'];
    unique_numbers(:,2) = unique_numbers(:,2) + base_num;
    
    % new numbers
    for istage = 1:numel(subj_corr_cell_cellid_predict{isubj_predict})
        for iu = 1:size(unique_numbers,1)
            subj_corr_cell_cellid_predict{isubj_predict}{istage}(subj_corr_cell_cellid_predict{isubj_predict}{istage}==unique_numbers(iu,1)) = unique_numbers(iu,2);
        end
    end
    
    % update base
    base_num = max(unique_numbers(:,2));
end

for isubj_unpredict = 1:length(unpredict_subjs)
    
    % unique cell numbers
    unique_numbers = unique(cell2mat(subj_corr_cell_cellid_unpredict{isubj_unpredict}(:)));
    unique_numbers = [unique_numbers (1:length(unique_numbers))'];
    unique_numbers(:,2) = unique_numbers(:,2) + base_num;
    
    % new numbers
    for istage = 1:numel(subj_corr_cell_cellid_unpredict{isubj_unpredict})
        for iu = 1:size(unique_numbers,1)
            subj_corr_cell_cellid_unpredict{isubj_unpredict}{istage}(subj_corr_cell_cellid_unpredict{isubj_unpredict}{istage}==unique_numbers(iu,1)) = unique_numbers(iu,2);
        end
    end
    
    % update base
    base_num = max(unique_numbers(:,2));
end

% days of interest
preprobe_session_numbers = 1:3:16;
firstday_session_numbers = 2:3:17;
lastday_session_numbers = 3:3:18;
postprobe_session_numbers = 4:3:19;

% preallcate subjects x probes
corr_predict_PreprobeProblem = cell(length(predict_subjs), length(preprobe_session_numbers));
    subject_idx_predict_PreprobeProblem = cell(size(corr_predict_PreprobeProblem));
    cell_id_predict_PreprobeProblem = cell(size(corr_predict_PreprobeProblem));
    stage_id_predict_PreprobeProblem = cell(size(corr_predict_PreprobeProblem));
    group_id_predict_PreprobeProblem = cell(size(corr_predict_PreprobeProblem));

corr_unpredict_PreprobeProblem = cell(length(unpredict_subjs), length(preprobe_session_numbers));
    subject_idx_unpredict_PreprobeProblem = cell(size(corr_unpredict_PreprobeProblem));
    cell_id_unpredict_PreprobeProblem = cell(size(corr_unpredict_PreprobeProblem));
    stage_id_unpredict_PreprobeProblem = cell(size(corr_unpredict_PreprobeProblem));
    group_id_unpredict_PreprobeProblem = cell(size(corr_unpredict_PreprobeProblem));
    
corr_predict_ProblemFinalprobe = cell(length(predict_subjs), length(preprobe_session_numbers));
    subject_idx_predict_ProblemFinalprobe = cell(size(corr_predict_ProblemFinalprobe));
    cell_id_predict_ProblemFinalprobe = cell(size(corr_predict_ProblemFinalprobe));
    stage_id_predict_ProblemFinalprobe = cell(size(corr_predict_ProblemFinalprobe));
    group_id_predict_ProblemFinalprobe = cell(size(corr_predict_ProblemFinalprobe));
    
corr_unpredict_ProblemFinalprobe = cell(length(unpredict_subjs), length(preprobe_session_numbers));
    subject_idx_unpredict_ProblemFinalprobe = cell(size(corr_unpredict_ProblemFinalprobe));
    cell_id_unpredict_ProblemFinalprobe = cell(size(corr_unpredict_ProblemFinalprobe));
    stage_id_unpredict_ProblemFinalprobe = cell(size(corr_unpredict_ProblemFinalprobe));
    group_id_unpredict_ProblemFinalprobe = cell(size(corr_unpredict_ProblemFinalprobe));
    
corr_predict_PreprobeFinalprobe = cell(length(predict_subjs), length(preprobe_session_numbers));
    subject_idx_predict_PreprobeFinalprobe = cell(size(corr_predict_PreprobeFinalprobe));
    cell_id_predict_PreprobeFinalprobe = cell(size(corr_predict_PreprobeFinalprobe));
    stage_id_predict_PreprobeFinalprobe = cell(size(corr_predict_PreprobeFinalprobe));
    group_id_predict_PreprobeFinalprobe = cell(size(corr_predict_PreprobeFinalprobe)); 
    
corr_unpredict_PreprobeFinalprobe = cell(length(unpredict_subjs), length(preprobe_session_numbers));
    subject_idx_unpredict_PreprobeFinalprobe = cell(size(corr_unpredict_PreprobeFinalprobe));
    cell_id_unpredict_PreprobeFinalprobe = cell(size(corr_unpredict_PreprobeFinalprobe));
    stage_id_unpredict_PreprobeFinalprobe = cell(size(corr_unpredict_PreprobeFinalprobe));
    group_id_unpredict_PreprobeFinalprobe = cell(size(corr_unpredict_PreprobeFinalprobe));
    
    
% iterate through sessions
for isesh = 1:length(preprobe_session_numbers)
    
    % adjacent probe sessions
    %
        % iterate through predict subjects
        for isubj_predict = 1:length(predict_subjs)
            % r values
            first_day = subj_corr_cell_predict{isubj_predict}{preprobe_session_numbers(isesh), firstday_session_numbers(isesh)};
            last_day = subj_corr_cell_predict{isubj_predict}{preprobe_session_numbers(isesh), lastday_session_numbers(isesh)};
            corr_predict_PreprobeProblem{isubj_predict, isesh} = [first_day; last_day];
            % subject number
            subject_idx_predict_PreprobeProblem{isubj_predict, isesh} = repmat(isubj_predict, size(corr_predict_PreprobeProblem{isubj_predict, isesh}));
            % cell number
            first_day = subj_corr_cell_cellid_predict{isubj_predict}{preprobe_session_numbers(isesh), firstday_session_numbers(isesh)};
            last_day = subj_corr_cell_cellid_predict{isubj_predict}{preprobe_session_numbers(isesh), lastday_session_numbers(isesh)};
            cell_id_predict_PreprobeProblem{isubj_predict, isesh} = [first_day; last_day];
            % stage number
            stage_id_predict_PreprobeProblem{isubj_predict, isesh} = repmat(isesh, size(corr_predict_PreprobeProblem{isubj_predict, isesh}));
            % group id
            group_id_predict_PreprobeProblem{isubj_predict, isesh} = zeros(size(corr_predict_PreprobeProblem{isubj_predict, isesh}));
        end
        
        % iterate through unpredict subjects
        for isubj_unpredict = 1:length(unpredict_subjs)
            % r values
            first_day = subj_corr_cell_unpredict{isubj_unpredict}{preprobe_session_numbers(isesh), firstday_session_numbers(isesh)};
            last_day = subj_corr_cell_unpredict{isubj_unpredict}{preprobe_session_numbers(isesh), lastday_session_numbers(isesh)};
            corr_unpredict_PreprobeProblem{isubj_unpredict, isesh} = [first_day; last_day];
            % subject number
            subject_idx_unpredict_PreprobeProblem{isubj_unpredict, isesh} = repmat(length(predict_subjs)+isubj_unpredict, size(corr_unpredict_PreprobeProblem{isubj_unpredict, isesh}));
            % cell number
            first_day = subj_corr_cell_cellid_unpredict{isubj_unpredict}{preprobe_session_numbers(isesh), firstday_session_numbers(isesh)};
            last_day = subj_corr_cell_cellid_unpredict{isubj_unpredict}{preprobe_session_numbers(isesh), lastday_session_numbers(isesh)};
            cell_id_unpredict_PreprobeProblem{isubj_unpredict, isesh} = [first_day; last_day];
            % stage number
            stage_id_unpredict_PreprobeProblem{isubj_unpredict, isesh} = repmat(isesh, size(corr_unpredict_PreprobeProblem{isubj_unpredict, isesh}));
            % group id
            group_id_unpredict_PreprobeProblem{isubj_unpredict, isesh} = ones(size(corr_unpredict_PreprobeProblem{isubj_unpredict, isesh}));
        end
        
    % compare all preprobes with last probe
    %
        % iterate through predict subjects
        for isubj_predict = 1:length(predict_subjs)
            % r values
            corr_predict_PreprobeFinalprobe{isubj_predict, isesh} = subj_corr_cell_predict{isubj_predict}{preprobe_session_numbers(isesh), postprobe_session_numbers(end)};
            % subject number
            subject_idx_predict_PreprobeFinalprobe{isubj_predict, isesh} = repmat(isubj_predict, size(corr_predict_PreprobeFinalprobe{isubj_predict, isesh}));
            % cell number
            cell_id_predict_PreprobeFinalprobe{isubj_predict, isesh} = subj_corr_cell_cellid_predict{isubj_predict}{preprobe_session_numbers(isesh), postprobe_session_numbers(end)};
            % stage number
            stage_id_predict_PreprobeFinalprobe{isubj_predict, isesh} = repmat(isesh, size(corr_predict_PreprobeFinalprobe{isubj_predict, isesh}));
            % group id
            group_id_predict_PreprobeFinalprobe{isubj_predict, isesh} = zeros(size(corr_predict_PreprobeFinalprobe{isubj_predict, isesh}));
        end
        
        % iterate through unpredict subjects
        for isubj_unpredict = 1:length(unpredict_subjs)
            % r values
            corr_unpredict_PreprobeFinalprobe{isubj_unpredict, isesh} = subj_corr_cell_unpredict{isubj_unpredict}{preprobe_session_numbers(isesh), postprobe_session_numbers(end)};
            % subject number
            subject_idx_unpredict_PreprobeFinalprobe{isubj_unpredict, isesh} = repmat(length(predict_subjs)+isubj_unpredict, size(corr_unpredict_PreprobeFinalprobe{isubj_unpredict, isesh}));
            % cell number
            cell_id_unpredict_PreprobeFinalprobe{isubj_unpredict, isesh} = subj_corr_cell_cellid_unpredict{isubj_unpredict}{preprobe_session_numbers(isesh), postprobe_session_numbers(end)};
            % stage number
            stage_id_unpredict_PreprobeFinalprobe{isubj_unpredict, isesh} = repmat(isesh, size(corr_unpredict_PreprobeFinalprobe{isubj_unpredict, isesh}));
            % group id
            group_id_unpredict_PreprobeFinalprobe{isubj_unpredict, isesh} = ones(size(corr_unpredict_PreprobeFinalprobe{isubj_unpredict, isesh}));
        end
        
    % compare all problems with last probe
    %
        % iterate through predict subjects
        for isubj_predict = 1:length(predict_subjs)
            % r values
            first_day = subj_corr_cell_predict{isubj_predict}{firstday_session_numbers(isesh), postprobe_session_numbers(end)};
            last_day = subj_corr_cell_predict{isubj_predict}{lastday_session_numbers(isesh), postprobe_session_numbers(end)};
            corr_predict_ProblemFinalprobe{isubj_predict, isesh} = [first_day; last_day];
            % subject number
            subject_idx_predict_ProblemFinalprobe{isubj_predict, isesh} = repmat(isubj_predict, size(corr_predict_ProblemFinalprobe{isubj_predict, isesh}));
            % cell number
            first_day = subj_corr_cell_cellid_predict{isubj_predict}{firstday_session_numbers(isesh), postprobe_session_numbers(end)};
            last_day = subj_corr_cell_cellid_predict{isubj_predict}{lastday_session_numbers(isesh), postprobe_session_numbers(end)};
            cell_id_predict_ProblemFinalprobe{isubj_predict, isesh} = [first_day; last_day];
            % stage number
            stage_id_predict_ProblemFinalprobe{isubj_predict, isesh} = repmat(isesh, size(corr_predict_ProblemFinalprobe{isubj_predict, isesh}));
            % group id
            group_id_predict_ProblemFinalprobe{isubj_predict, isesh} = zeros(size(corr_predict_ProblemFinalprobe{isubj_predict, isesh}));
        end
        
        % iterate through unpredict subjects
        for isubj_unpredict = 1:length(unpredict_subjs)
            % r values
            first_day = subj_corr_cell_unpredict{isubj_unpredict}{firstday_session_numbers(isesh), postprobe_session_numbers(end)};
            last_day = subj_corr_cell_unpredict{isubj_unpredict}{lastday_session_numbers(isesh), postprobe_session_numbers(end)};
            corr_unpredict_ProblemFinalprobe{isubj_unpredict, isesh} = [first_day; last_day];
            % subject number
            subject_idx_unpredict_ProblemFinalprobe{isubj_unpredict, isesh} = repmat(length(predict_subjs)+isubj_unpredict, size(corr_unpredict_ProblemFinalprobe{isubj_unpredict, isesh}));
            % cell number
            first_day = subj_corr_cell_cellid_unpredict{isubj_unpredict}{firstday_session_numbers(isesh), postprobe_session_numbers(end)};
            last_day = subj_corr_cell_cellid_unpredict{isubj_unpredict}{lastday_session_numbers(isesh), postprobe_session_numbers(end)};
            cell_id_unpredict_ProblemFinalprobe{isubj_unpredict, isesh} = [first_day; last_day];
            % stage number
            stage_id_unpredict_ProblemFinalprobe{isubj_unpredict, isesh} = repmat(isesh, size(corr_unpredict_ProblemFinalprobe{isubj_unpredict, isesh}));
            % group id
            group_id_unpredict_ProblemFinalprobe{isubj_unpredict, isesh} = ones(size(corr_unpredict_ProblemFinalprobe{isubj_unpredict, isesh}));            
        end
end

% combine corrlation information for later mixed models
datamtx_PreprobeProblem = [cell2mat(corr_predict_PreprobeProblem(:)); cell2mat(corr_unpredict_PreprobeProblem(:))];
    subjid_PreprobeProblem = [cell2mat(subject_idx_predict_PreprobeProblem(:)); cell2mat(subject_idx_unpredict_PreprobeProblem(:))];
    cellid_PreprobeProblem = [cell2mat(cell_id_predict_PreprobeProblem(:)); cell2mat(cell_id_unpredict_PreprobeProblem(:))];
    stageid_PreprobeProblem = [cell2mat(stage_id_predict_PreprobeProblem(:)); cell2mat(stage_id_unpredict_PreprobeProblem(:))];
    groupid_PreprobeProblem = [cell2mat(group_id_predict_PreprobeProblem(:)); cell2mat(group_id_unpredict_PreprobeProblem(:))];
datamtx_PreprobeFinalprobe = [cell2mat(corr_predict_PreprobeFinalprobe(:)); cell2mat(corr_unpredict_PreprobeFinalprobe(:))];
    subjid_PreprobeFinalprobe = [cell2mat(subject_idx_predict_PreprobeFinalprobe(:)); cell2mat(subject_idx_unpredict_PreprobeFinalprobe(:))];
    cellid_PreprobeFinalprobe = [cell2mat(cell_id_predict_PreprobeFinalprobe(:)); cell2mat(cell_id_unpredict_PreprobeFinalprobe(:))];
    stageid_PreprobeFinalprobe = [cell2mat(stage_id_predict_PreprobeFinalprobe(:)); cell2mat(stage_id_unpredict_PreprobeFinalprobe(:))];
    groupid_PreprobeFinalprobe = [cell2mat(group_id_predict_PreprobeFinalprobe(:)); cell2mat(group_id_unpredict_PreprobeFinalprobe(:))];    
datamtx_ProblemFinalprobe = [cell2mat(corr_predict_ProblemFinalprobe(:)); cell2mat(corr_unpredict_ProblemFinalprobe(:))];
    subjid_ProblemFinalprobe = [cell2mat(subject_idx_predict_ProblemFinalprobe(:)); cell2mat(subject_idx_unpredict_ProblemFinalprobe(:))];
    cellid_ProblemFinalprobe = [cell2mat(cell_id_predict_ProblemFinalprobe(:)); cell2mat(cell_id_unpredict_ProblemFinalprobe(:))];
    stageid_ProblemFinalprobe = [cell2mat(stage_id_predict_ProblemFinalprobe(:)); cell2mat(stage_id_unpredict_ProblemFinalprobe(:))];
    groupid_ProblemFinalprobe = [cell2mat(group_id_predict_ProblemFinalprobe(:)); cell2mat(group_id_unpredict_ProblemFinalprobe(:))];    

    
    
    
    
    
% errorbar plots
%

    % mixed model   
    model_str = 'PopCorr~Session*Group+(1|Subject)+(1|Neuron)';
    tbl = table(datamtx_PreprobeProblem, stageid_PreprobeProblem, groupid_PreprobeProblem, subjid_PreprobeProblem, cellid_PreprobeProblem, 'VariableNames', {'PopCorr', 'Session', 'Group', 'Subject', 'Neuron'});
    tbl.Subject = categorical(tbl.Subject);
    tbl.Neuron = categorical(tbl.Neuron);
    lme_PreprobeProblem = fitlme(tbl, model_str)
    
            % PreprobeProblem plot 
            plot_cell_predict = cell(1, length(unique(stageid_PreprobeProblem)));
            plot_cell_unpredict = cell(1, length(unique(stageid_PreprobeProblem)));
                plot_cell_predict_subjmeans = nan(length(predict_subjs), length(unique(stageid_PreprobeProblem)));
                plot_cell_unpredict_subjmeans = nan(length(predict_subjs), length(unique(stageid_PreprobeProblem)));
            for istage = 1:length(unique(stageid_PreprobeProblem))
                plot_cell_predict{istage} = datamtx_PreprobeProblem(stageid_PreprobeProblem==istage & groupid_PreprobeProblem==0);
                plot_cell_unpredict{istage} = datamtx_PreprobeProblem(stageid_PreprobeProblem==istage & groupid_PreprobeProblem==1);
                for isubj = 1:length(predict_subjs)
                    plot_cell_predict_subjmeans(isubj, istage) = mean(datamtx_PreprobeProblem(stageid_PreprobeProblem==istage & groupid_PreprobeProblem==0 & subjid_PreprobeProblem==isubj));
                end
                for isubj = length(predict_subjs)+1:length(predict_subjs)+length(unpredict_subjs)
                    plot_cell_unpredict_subjmeans(isubj, istage) = mean(datamtx_PreprobeProblem(stageid_PreprobeProblem==istage & groupid_PreprobeProblem==1 & subjid_PreprobeProblem==isubj));
                end
            end
            figure; hold on
            errorbar_plot(plot_cell_unpredict, 0, [], light_blue_color, blue_color)
            errorbar_plot(plot_cell_predict, 0, [], light_green_color, green_color)
            tic_vect = [-.999 -.99 -.9 -.8 -.6 -.3 0 .3 .6 .8 .9 .99 .999]; ylim(atanh([-.99 .99])); yticks(atanh(tic_vect)); yticklabels(tic_vect)
            xticks(1:6); xticklabels(0:6);
            title('Preprobe-Problem')
            %figure; hold on;
            %errorbar_mtx(plot_cell_predict_subjmeans)
            %errorbar_mtx(plot_cell_unpredict_subjmeans)
            
    
    model_str = 'PopCorr~Session*Group+(1|Subject)+(1|Neuron)';
    tbl = table(datamtx_PreprobeFinalprobe, stageid_PreprobeFinalprobe, groupid_PreprobeFinalprobe, subjid_PreprobeFinalprobe, cellid_PreprobeFinalprobe, 'VariableNames', {'PopCorr', 'Session', 'Group', 'Subject', 'Neuron'});
    tbl.Subject = categorical(tbl.Subject);
    tbl.Neuron = categorical(tbl.Neuron);
    lme_PreprobeFinalprobe = fitlme(tbl, model_str)
    
            % PreprobeFinalprobe plot 
            plot_cell_predict = cell(1, length(unique(stageid_PreprobeFinalprobe)));
            plot_cell_unpredict = cell(1, length(unique(stageid_PreprobeFinalprobe)));
                plot_cell_predict_subjmeans = nan(length(predict_subjs), length(unique(stageid_PreprobeFinalprobe)));
                plot_cell_unpredict_subjmeans = nan(length(predict_subjs), length(unique(stageid_PreprobeFinalprobe)));
            for istage = 1:length(unique(stageid_PreprobeFinalprobe))
                plot_cell_predict{istage} = datamtx_PreprobeFinalprobe(stageid_PreprobeFinalprobe==istage & groupid_PreprobeFinalprobe==0);
                plot_cell_unpredict{istage} = datamtx_PreprobeFinalprobe(stageid_PreprobeFinalprobe==istage & groupid_PreprobeFinalprobe==1);
                for isubj = 1:length(predict_subjs)
                    plot_cell_predict_subjmeans(isubj, istage) = mean(datamtx_PreprobeFinalprobe(stageid_PreprobeFinalprobe==istage & groupid_PreprobeFinalprobe==0 & subjid_PreprobeFinalprobe==isubj));
                end
                for isubj = length(predict_subjs)+1:length(predict_subjs)+length(unpredict_subjs)
                    plot_cell_unpredict_subjmeans(isubj, istage) = mean(datamtx_PreprobeFinalprobe(stageid_PreprobeFinalprobe==istage & groupid_PreprobeFinalprobe==1 & subjid_PreprobeFinalprobe==isubj));
                end
            end
            figure; hold on
            errorbar_plot(plot_cell_unpredict, 0, [], light_blue_color, blue_color)
            errorbar_plot(plot_cell_predict, 0, [], light_green_color, green_color)
            tic_vect = [-.999 -.99 -.9 -.8 -.6 -.3 0 .3 .6 .8 .9 .99 .999]; ylim(atanh([-.99 .99])); yticks(atanh(tic_vect)); yticklabels(tic_vect)
            xticks(1:6); xticklabels(0:6);
            title('Preprobe-FinalProbe')
            %figure; hold on;
            %errorbar_mtx(plot_cell_predict_subjmeans)
            %errorbar_mtx(plot_cell_unpredict_subjmeans)

    model_str = 'PopCorr~Session*Group+(1|Subject)+(1|Neuron)';
    tbl = table(datamtx_ProblemFinalprobe, stageid_ProblemFinalprobe, groupid_ProblemFinalprobe, subjid_ProblemFinalprobe, cellid_ProblemFinalprobe, 'VariableNames', {'PopCorr', 'Session', 'Group', 'Subject', 'Neuron'});
    tbl.Subject = categorical(tbl.Subject);
    tbl.Neuron = categorical(tbl.Neuron);
    lme_ProblemFinalprobe = fitlme(tbl, model_str)
    
            % ProblemFinalprobe plot 
            plot_cell_predict = cell(1, length(unique(stageid_ProblemFinalprobe)));
            plot_cell_unpredict = cell(1, length(unique(stageid_ProblemFinalprobe)));
                plot_cell_predict_subjmeans = nan(length(predict_subjs), length(unique(stageid_ProblemFinalprobe)));
                plot_cell_unpredict_subjmeans = nan(length(predict_subjs), length(unique(stageid_ProblemFinalprobe)));
            for istage = 1:length(unique(stageid_ProblemFinalprobe))
                plot_cell_predict{istage} = datamtx_ProblemFinalprobe(stageid_ProblemFinalprobe==istage & groupid_ProblemFinalprobe==0);
                plot_cell_unpredict{istage} = datamtx_ProblemFinalprobe(stageid_ProblemFinalprobe==istage & groupid_ProblemFinalprobe==1);
                for isubj = 1:length(predict_subjs)
                    plot_cell_predict_subjmeans(isubj, istage) = mean(datamtx_ProblemFinalprobe(stageid_ProblemFinalprobe==istage & groupid_ProblemFinalprobe==0 & subjid_ProblemFinalprobe==isubj));
                end
                for isubj = length(predict_subjs)+1:length(predict_subjs)+length(unpredict_subjs)
                    plot_cell_unpredict_subjmeans(isubj, istage) = mean(datamtx_ProblemFinalprobe(stageid_ProblemFinalprobe==istage & groupid_ProblemFinalprobe==1 & subjid_ProblemFinalprobe==isubj));
                end
            end
            figure; hold on
            errorbar_plot(plot_cell_unpredict, 0, [], light_blue_color, blue_color)
            errorbar_plot(plot_cell_predict, 0, [], light_green_color, green_color)
            tic_vect = [-.999 -.99 -.9 -.8 -.6 -.3 0 .3 .6 .8 .9 .99 .999]; ylim(atanh([-.99 .99])); yticks(atanh(tic_vect)); yticklabels(tic_vect)
            xticks(1:6); xticklabels(0:6);
            title('Problem-FinalProbe')
            %figure; hold on;
            %errorbar_mtx(plot_cell_predict_subjmeans)
            %errorbar_mtx(plot_cell_unpredict_subjmeans)

%}


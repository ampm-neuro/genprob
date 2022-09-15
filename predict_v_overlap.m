% prediction scores vs overlap

training_group_folders = [{'train_mevar_imaging_hpc_fin'} {'train_hivar_imaging_hpc_fin'}];

for itgf = 1:length(training_group_folders)
    training_group_folder = training_group_folders{itgf};

    %% unique subject files for this training group
    training_group_folder_path = ['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\' training_group_folder];
    subject_folders = get_folder_paths_all(training_group_folder_path);
    unq_subjs = cell(size(subject_folders));
    for isubj = 1:size(subject_folders,1)
        unq_subjs{isubj} = find_subj_id(subject_folders{isubj});
    end


    %% Problem tones for this training group

    % all tones
    load('unqfrq41', 'unqfrq41')

    % tones in each problem (rich, poor)
    if contains(training_group_folder, 'mevar')
        all_prob_tones = rich_bounds_prob('mevar', 0);
    elseif contains(training_group_folder, 'hivar')
        all_prob_tones = rich_bounds_prob('hivar', 0);
    else
        error('input error')
    end

    % number of problems
    num_problems = size(all_prob_tones,1);

    % tones numbers (1:41) in each problem (rich, poor)
    all_prob_tone_nums = nan(size(all_prob_tones));
    for iprob = 1:num_problems
        for itone = 1:size(all_prob_tones,2)
            all_prob_tone_nums(iprob, itone) = find(ismember(unqfrq41, all_prob_tones(iprob, itone)));
        end
    end


    %% all probe wait times

    % wait times from every probe for every subject (probe,tone,subj)
    [probe_wait_times, unq_subjects_probe_all, subj_path_cell] = ALL_probe_wait_times(training_group_folder, 0); % observed values only

    % number of probes
    num_probes = size(probe_wait_times,1);
    num_tones = size(probe_wait_times,2);
    num_subjects = size(probe_wait_times,3);



    %% Fit model to each probe session (wait times)
    %
    % preallcoate 
    num_coefs = 4;
    coefficients = nan(num_probes, num_coefs, num_subjects);
    subj_curves = nan(num_probes, length(1:0.001:num_tones), num_subjects);

    % iterate through all probes
    for iprobe = 1:num_probes

        % open curve of all subject curves for this probe
        figure; hold on

        % iterate through subjects
        for isubj = 1:num_subjects

            % isolate probe waits
            wait_times = probe_wait_times(iprobe, :, isubj);

            % compute coefficients
            try
                mean_wait_times = nanmean(probe_wait_times(iprobe,:,isubj));
                %[~, coefEsts, modelFun] = ampm_normal_logistic_fit(1:num_tones, wait_times, [mean_wait_times 21 0 1]);
                [~, coefEsts, modelFun] = ampm_normal_logistic_fit_algo(1:num_tones, wait_times, [mean_wait_times 20 0 1]);

            % if modelling fails, set to mean
            catch 

                coefEsts = [nanmean(probe_wait_times(iprobe,:,isubj)) zeros(1,num_coefs-1)];
            end

            % load coefficients
            coefficients(iprobe,:,isubj) = coefEsts;

            % santity check coefs
            %acceptable_coef_bounds = [-5 45; -5 41; -20 20; -20 20];
            acceptable_coef_bounds = [-30 60; 1 41; -30 30; -30 30];
            if any(...
                    coefficients(iprobe,1,isubj) < acceptable_coef_bounds(1,1) | coefficients(iprobe,1,isubj) > acceptable_coef_bounds(1,2) | isnan(coefficients(iprobe,1,isubj)) | ...
                    coefficients(iprobe,2,isubj) < acceptable_coef_bounds(2,1) | coefficients(iprobe,2,isubj) > acceptable_coef_bounds(2,2) | isnan(coefficients(iprobe,2,isubj)) |...
                    coefficients(iprobe,3,isubj) < acceptable_coef_bounds(3,1) | coefficients(iprobe,3,isubj) > acceptable_coef_bounds(3,2) | isnan(coefficients(iprobe,3,isubj)) |...
                    coefficients(iprobe,4,isubj) < acceptable_coef_bounds(4,1) | coefficients(iprobe,4,isubj) > acceptable_coef_bounds(4,2) | isnan(coefficients(iprobe,4,isubj)) ...
                    )
                coefficients(iprobe,:,isubj) = [nanmean(probe_wait_times(iprobe,:,isubj)) zeros(1,num_coefs-1)];
            end

            % deal with missing data
            if isnan(coefficients(iprobe,1,isubj))
                coefficients(iprobe,:,isubj) = nan;
            end

        end

    end


    %% Use coefficients to estimate discrimination of every problem tone pair during every probe test

    %
    % preallocate (probe, problem, subj)
    rich_waits_problem = nan(num_probes, num_problems, num_subjects);
    poor_waits_problem = nan(num_probes, num_problems, num_subjects);
    diff_waits_problem = nan(num_probes, num_problems, num_subjects);

    % iterate through all probes
    for isubj = 1:size(probe_wait_times,3)
        for iproblem = 1:num_problems
            for iprobe = 1:num_probes

                % tone numbers
                rich_tone_num = all_prob_tone_nums(iproblem, 1);
                poor_tone_num = all_prob_tone_nums(iproblem, 2);

                % tone wait times (during probe) FULL MODEL
                rich_waits_problem(iprobe, iproblem, isubj) = modelFun(coefficients(iprobe,:,isubj), rich_tone_num);
                poor_waits_problem(iprobe, iproblem, isubj) = modelFun(coefficients(iprobe,:,isubj), poor_tone_num);

                % difference between rich and poor tone wait times (during probe)
                diff_waits_problem(iprobe, iproblem, isubj)...
                    = rich_waits_problem(iprobe, iproblem, isubj) - poor_waits_problem(iprobe, iproblem, isubj);
            end
        end
    end



    %% discrimination of previous problem during each probe
    diff_waits_problem_preprobe = nan(num_subjects, num_problems);
    diff_waits_problem_postprobe = nan(num_subjects, num_problems);
    diff_waits_problem_2postprobe = nan(num_subjects, num_problems-1);
    all_other_past_probes = nan(num_subjects, num_problems);
    all_future_probes = nan(num_subjects, num_problems);
    for isubj = 1:num_subjects

        % for current subject
        subj_dwp_preprobe = diff_waits_problem(1:6, :, isubj);
        subj_dwp_postprobe = diff_waits_problem(2:7, :, isubj);
        subj_dwp_2postprobe = diff_waits_problem(3:7, :, isubj);

        % extract probe and problem matches
        diff_waits_problem_preprobe(isubj,:) = subj_dwp_preprobe(logical(eye(size(subj_dwp_preprobe))));
        diff_waits_problem_postprobe(isubj,:) = subj_dwp_postprobe(logical(eye(size(subj_dwp_postprobe))));
        diff_waits_problem_2postprobe(isubj,:) = subj_dwp_2postprobe(logical(eye(size(subj_dwp_2postprobe))));

        % extract additional probe and problem matches
        for iprobe = 3:7
            all_other_past_probes(isubj, iprobe-1) = nanmean(diff_waits_problem(iprobe, 1:iprobe-2, isubj));
        end
        for iprobe = 1:6
            all_future_probes(isubj, iprobe) = nanmean(diff_waits_problem(iprobe, iprobe:end, isubj));
        end


    end

    %% all future problems prediction score
    future_dwp = diff_waits_problem;
    future_dwp = future_dwp(1:6,:,:);
    for isubj = 1:size(diff_waits_problem,3)
        for iprobe = 1:6
            for iproblem = 1:6
                if iprobe>iproblem
                    future_dwp(iprobe, iproblem, isubj) = nan;
                end
            end
        end
    end



    %% overlap values
    preprobe_idx = 1:3:16;
    postprobe_idx = 4:3:19;
    problem_firstday_idx = preprobe_idx+1;
    problem_secondday_idx = preprobe_idx+2;
    overlap_mtx = overlapping_population_average(unq_subjs);

    % next problem
    %
    overlap_mtx_idx = sub2ind([20 20], problem_firstday_idx, preprobe_idx);
    overlap_next_problem_vals = [];
    for isubj = 1:size(overlap_mtx,3)
        overlap_mtx_local = overlap_mtx(:,:,isubj);
        subj_overlap_vals = overlap_mtx_local(overlap_mtx_idx);
        overlap_next_problem_vals = [overlap_next_problem_vals subj_overlap_vals'];
    end
    overlap_next_problem_vals = overlap_next_problem_vals';
    %}

    % all future problems
    %{
    overlap_mtx_idx = sub2ind([20 20], preprobe_idx+1, preprobe_idx);
    overlap_future_problem_vals = nan(size(future_dwp));
    for isubj = 1:size(overlap_mtx,3)
        for iprobe = 1:6
            for iproblem = 1:6
                if iproblem<iprobe
                    continue
                end
                overlap_future_problem_vals(iprobe, iproblem, isubj) =...
                    overlap_mtx(problem_firstday_idx(iproblem), preprobe_idx(iprobe), isubj);
            end
        end
    end
    %}


    % next probe
    %
    overlap_mtx_idx = sub2ind([20 20], postprobe_idx, preprobe_idx);
    overlap_next_probe_vals = [];
    for isubj = 1:size(overlap_mtx,3)
        overlap_mtx_local = overlap_mtx(:,:,isubj);
        subj_overlap_vals = overlap_mtx_local(overlap_mtx_idx);
        overlap_next_probe_vals = [overlap_next_probe_vals subj_overlap_vals'];
    end
    overlap_next_probe_vals = overlap_next_probe_vals';
    %}



    %% quality of prediction related to overlap

        accuracy_of_prediction = diff_waits_problem_preprobe;
        overlap_problems = overlap_next_problem_vals;
        overlap_probes = overlap_next_probe_vals;
        %accuracy_of_prediction = future_dwp;
        %overlap = overlap_future_problem_vals;

        figure; 
        subplot(1,2,1); hold on;
        [r, p] = fit_line(accuracy_of_prediction(:), overlap_problems(:));
        plot(xlim, [0 0], 'k--'); plot([0 0], ylim, 'k--')
        title NextProblem
        subplot(1,2,2); hold on;
        [r, p] = fit_line(accuracy_of_prediction(:), overlap_probes(:));
        plot(xlim, [0 0], 'k--'); plot([0 0], ylim, 'k--')
        title NextProbe


    % tones in each problem (rich, poor)
    if contains(training_group_folder, 'mevar')
        accuracy_of_prediction_mevar = accuracy_of_prediction;
        overlap_problems_mevar = overlap_problems;
        overlap_probes_mevar = overlap_probes;
    elseif contains(training_group_folder, 'hivar')
        accuracy_of_prediction_hivar = accuracy_of_prediction;
        overlap_problems_hivar = overlap_problems;
        overlap_probes_hivar = overlap_probes;
    end


    

end

%% combine mevar and hivar
try
    accuracy_of_prediction_all = cat(3, accuracy_of_prediction_mevar, accuracy_of_prediction_hivar);
    overlap_problems_all = cat(3, overlap_problems_mevar, overlap_problems_hivar);
    overlap_probes_all = cat(3, overlap_probes_mevar, overlap_probes_hivar);
catch
    accuracy_of_prediction_all = [accuracy_of_prediction_mevar; accuracy_of_prediction_hivar];
    overlap_problems_all = [overlap_problems_mevar; overlap_problems_hivar];
    overlap_probes_all = [overlap_probes_mevar; overlap_probes_hivar];
end


%% overall plot
figure; 

subplot(1,2,1); hold on;
plot(accuracy_of_prediction_mevar, overlap_problems_mevar, 'go')
plot(accuracy_of_prediction_hivar, overlap_problems_hivar, 'bo')
[r, p] = fit_line(accuracy_of_prediction_all(:), overlap_problems_all(:), 0);
title(['Probe-to-problem; r=' num2str(r) '; p=' num2str(p)])
axis([-20 20 -.05 .6]); axis square
plot(xlim, [0 0], 'k--'); plot([0 0], ylim, 'k--')
xlabel('Prediction accuracy (s)')
ylabel('Proportion overlap')

subplot(1,2,2); hold on;
plot(accuracy_of_prediction_mevar, overlap_probes_mevar, 'go')
plot(accuracy_of_prediction_hivar, overlap_probes_hivar, 'bo')
[r, p] = fit_line(accuracy_of_prediction_all(:), overlap_probes_all(:), 0);
title(['Probe-to-probe; r=' num2str(r) '; p=' num2str(p)])
axis([-20 20 -.05 .6]); axis square
plot(xlim, [0 0], 'k--'); plot([0 0], ylim, 'k--')
xlabel('Prediction accuracy (s)')
ylabel('Proportion overlap')

    





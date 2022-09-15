function [logist_coef_cell_corrected, norm_coef_cell] = ALL_model_probes(training_group, probe_nums, plot_OnOFF)
% Used to compare training problem behavior with probe behavior in
% every which way. produces many plots. additional ones found near the
% bottom


%% all probe wait times

% smooth wait times from every probe for every subject (probe,tone,subj)
[probe_wait_times] = ALL_probe_wait_times(training_group, 0);


% only input prob nums
probe_wait_times = probe_wait_times(probe_nums,:,:);

% number of probes
num_probes = size(probe_wait_times,1);
num_tones = size(probe_wait_times,2);
num_subjects = size(probe_wait_times,3);

% zscored probe wait times
probe_wait_times_z = nan(num_probes, num_tones, num_subjects);
for iprobe = 1:num_probes
    for isubj = 1:num_subjects
        probe_wait_times_z(iprobe, :, isubj) = zscore_mtx(probe_wait_times(iprobe, :, isubj)')';
    end
end



%% Fit model to each probe session (wait times)
%
% preallcoate 
num_coefs = 4;
coefficients = nan(num_probes, num_coefs, num_subjects);
subj_curves = nan(num_probes, length(1:0.001:num_tones), num_subjects);

% iterate through all probes
for iprobe = 1:num_probes

    % iterate through subjects
    for isubj = 1:num_subjects
        
        % isolate probe waits
        wait_times = probe_wait_times(iprobe, :, isubj);
        if sum(~isnan(wait_times))<10
            continue
        end
        
        % compute coefficients
        %try
            mean_wait_times = nanmean(probe_wait_times(iprobe,:,isubj));
            [~, coefEsts, modelFun] = ampm_normal_logistic_fit_algo(1:num_tones, wait_times, [mean_wait_times 20 0 1]);
            
         %if modelling fails, set to mean
        %catch 
        %    coefEsts = [nanmean(probe_wait_times(iprobe,:,isubj)) zeros(1,num_coefs-1)];
       % end
        
        % load coefficients
        coefficients(iprobe,:,isubj) = coefEsts;
        
        % santity check coefs
        if any(...
                coefficients(iprobe,1,isubj) > 45 | coefficients(iprobe,1,isubj) < -5 | isnan(coefficients(iprobe,1,isubj)) | ...
                coefficients(iprobe,2,isubj) > 41 | coefficients(iprobe,2,isubj) < -5 | isnan(coefficients(iprobe,1,isubj)) |...
                coefficients(iprobe,3,isubj) > 20 | coefficients(iprobe,3,isubj) < -20 | isnan(coefficients(iprobe,1,isubj)) |...
                coefficients(iprobe,4,isubj) > 20 | coefficients(iprobe,4,isubj) < -20 | isnan(coefficients(iprobe,1,isubj)) ...
                )
            coefficients(iprobe,:,isubj) = [nanmean(probe_wait_times(iprobe,:,isubj)) zeros(1,num_coefs-1)];
        end
        
        % deal with missing data
        if isnan(coefficients(iprobe,1,isubj))
            coefficients(iprobe,:,isubj) = nan;
        end
        
        
        if any(~isnan(wait_times))
            
            %model
            [subj_curves(iprobe,:,isubj)] = modelFun(coefEsts, (1:0.001:num_tones)');
            
            % plots
            if plot_OnOFF == 1

                figure; hold on
                plot(1:num_tones, wait_times, 'o')
                plot(1:0.001:num_tones, subj_curves(iprobe,:,isubj), '-', 'color', 0.7.*[1 1 1], 'linewidth', 2)

                    % overlay log and norm components
                    coefEsts_log = coefEsts; coefEsts_log(4) = realmin;
                    plot(1:0.001:num_tones, modelFun(coefEsts_log, (1:0.001:num_tones)'), '--', 'color', 'r', 'linewidth', 1)
                    coefEsts_norm = coefEsts; coefEsts_norm(1) = coefEsts_norm(1) + coefEsts_norm(3)/2; coefEsts_norm(3) = realmin;
                    plot(1:0.001:num_tones, modelFun(coefEsts_norm, (1:0.001:num_tones)'), '--', 'color', 'b', 'linewidth', 1)

                title(['Subj ' num2str(isubj) '; Probe ' num2str(probe_nums(iprobe)) '; Norm = ' num2str(coefficients(iprobe,4,isubj)) ', Logistic = ' num2str(coefficients(iprobe,3,isubj))])
                xlim([0.5 41.5])
                ylim([0 60])
                ylabel('Wait times (s)')
                xlabel('Tone number')
            end
        end
    end
    
    % plot mean and se of all subject curves
    %{
    mean_fit_curve = nanmean(subj_curves(iprobe,:,:),3);
    se_fit_curves = nanstd(subj_curves(iprobe,:,:),[],3)./sqrt(sum(~isnan(subj_curves(iprobe,1,:)),3));
    plot(1:0.001:num_tones, mean_fit_curve, 'k-', 'linewidth', 4)
    plot(1:0.001:num_tones, mean_fit_curve-se_fit_curves, 'k-', 'linewidth', 2)
    plot(1:0.001:num_tones, mean_fit_curve+se_fit_curves, 'k-', 'linewidth', 2)
    set(gca,'TickLength',[0, 0]); box off;
    
    title(['Probe ' num2str(iprobe) ': ])
    %}
    
end



%% Plot model coefficients

% Logistic curve mulitplier
logist_ceof = squeeze(coefficients(:,3,:))';
logist_ceof_corrected = squeeze(coefficients(:,3,:))';
logist_coef_cell = cell(1,9);
logist_coef_cell_corrected = cell(1,9);
for iprobe = 1:num_probes
    if ismember(iprobe, [2 4 6])
        logist_coef_cell_corrected{iprobe} = -logist_ceof(:,iprobe);
        logist_ceof_corrected(:,iprobe) = -logist_ceof(:,iprobe);
    else
        logist_coef_cell_corrected{iprobe} = logist_ceof(:,iprobe);
    end
    logist_coef_cell{iprobe} = logist_ceof(:,iprobe);
end
%
figure; hold on;
%errorbar_plot(logist_coef_cell_corrected(1), 0, 1);
%errorbar_plot(logist_coef_cell_corrected(2:7), 1, 2:7);
%errorbar_plot(logist_coef_cell_corrected(8), 0, 8);
%errorbar_plot(logist_coef_cell_corrected(9), 0, 9);
errorbar_plot(logist_coef_cell(1), 0, 1);
errorbar_plot(logist_coef_cell(2:7), 1, 2:7);
errorbar_plot(logist_coef_cell(8), 0, 8);
errorbar_plot(logist_coef_cell(9), 0, 9);
xlim([0.5 num_probes+0.5])
plot(xlim, [1 1].*0, 'k--')
set(gca,'TickLength',[0, 0]); box off;
xlabel('Probe number')
ylabel('Wait times (s)')
title('Logistic curve height')
%}

% Normal curve mulitplier
norm_coef = squeeze(coefficients(:,4,:))';
norm_coef_cell = cell(1,8);
for iprobe = 1:num_probes
    norm_coef_cell{iprobe} = norm_coef(:,iprobe);
end
%
figure; hold on;
errorbar_plot(norm_coef_cell(1), 0, 1);
errorbar_plot(norm_coef_cell(2:7), 1, 2:7);
errorbar_plot(norm_coef_cell(8), 0, 8);
errorbar_plot(norm_coef_cell(9), 0, 9);
xlim([0.5 num_probes+0.5])
plot(xlim, [1 1].*0, 'k--')
set(gca,'TickLength',[0, 0]); box off;
xlabel('Probe number')
ylabel('Wait times (s)')
title('Normal curve height')
%}

% Intercept of curve
intercept_coef = squeeze(coefficients(:,1,:))';
intercept_coef_cell = cell(1,9);
for iprobe = 1:num_probes
    intercept_coef_cell{iprobe} = intercept_coef(:,iprobe);
end
%{
figure; hold on;
errorbar_plot(intercept_coef_cell(1), 0, 1);
errorbar_plot(intercept_coef_cell(2:7), 1, 2:7);
errorbar_plot(intercept_coef_cell(8), 0, 8);
errorbar_plot(intercept_coef_cell(9), 0, 9);
xlim([0.5 num_probes+0.5])
set(gca,'TickLength',[0, 0]); box off;
xlabel('Probe number')
ylabel('Wait times (s)')
title('Intercept of fit curve')
%}

% Tonal center of curve
center_coef = squeeze(coefficients(:,3,:))';
center_coef_cell = cell(1,9);
for iprobe = 1:num_probes
    center_coef_cell{iprobe} = center_coef(:,iprobe);
end
%{
figure; hold on;
errorbar_plot(center_coef_cell(1), 0, 1);
errorbar_plot(center_coef_cell(2:7), 1, 2:7);
errorbar_plot(center_coef_cell(8), 0, 8);
errorbar_plot(center_coef_cell(9), 0, 9);
xlim([0.5 num_probes+0.5])
set(gca,'TickLength',[0, 0]); box off;
xlabel('Probe number')
ylabel('Curve centers (tone number)')
title('Center of fit curve')
%}








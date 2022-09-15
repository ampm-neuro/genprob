


%% get preprobe and first problem wait times
[train_mean_waits_mevar, probe_model_waits_mevar] = probe_predict_lastproblem('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_imaging_hpc');
[train_mean_waits_hivar, probe_model_waits_hivar] = probe_predict_lastproblem('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_hivar_imaging_hpc');


% prediction strength (rich - poor)
%
% training r-p - preprobe r-p ; therefore, high values mean prediction
% improved, low values mean prediction got worse, extreme values mean
% prediction changed in general
%
subj_mevar_predstr = nan(6,length(train_mean_waits_mevar));
for isubj_mevar = 1:length(train_mean_waits_mevar)
    subj_mevar_predstr(:,isubj_mevar) = (train_mean_waits_mevar{isubj_mevar}(:,1)-train_mean_waits_mevar{isubj_mevar}(:,2)) - (probe_model_waits_mevar{isubj_mevar}(:,1)-probe_model_waits_mevar{isubj_mevar}(:,2));
end
subj_hivar_predstr = nan(6,length(train_mean_waits_hivar));
for isubj_hivar = 1:length(train_mean_waits_hivar)
    subj_hivar_predstr(:,isubj_hivar) = (train_mean_waits_hivar{isubj_hivar}(:,1)-train_mean_waits_hivar{isubj_hivar}(:,2)) - (probe_model_waits_hivar{isubj_hivar}(:,1)-probe_model_waits_hivar{isubj_hivar}(:,2));
end

subj_mevar_predstr %(session, subj)
subj_hivar_predstr %(session, subj)



%% pre to post-probe coefficient changes
% get normal and log coefs for every subject and probe
[coefEsts_out_normal_mevar, coefEsts_out_log_mevar] = probe_coefs('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_imaging_hpc');
coefEsts_out_normal_mevar_delta = coefEsts_out_normal_mevar(2:end,:) - coefEsts_out_normal_mevar(1:end-1,:);
coefEsts_out_log_mevar_delta = coefEsts_out_log_mevar(2:end,:) - coefEsts_out_log_mevar(1:end-1,:);
[coefEsts_out_normal_hivar, coefEsts_out_log_hivar] = probe_coefs('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_hivar_imaging_hpc');
coefEsts_out_normal_hivar_delta = coefEsts_out_normal_hivar(2:end,:) - coefEsts_out_normal_hivar(1:end-1,:);
coefEsts_out_log_hivar_delta = coefEsts_out_log_hivar(2:end,:) - coefEsts_out_log_hivar(1:end-1,:);

% CORRECT LOG FOR DIRECTIONALITY
coefEsts_out_log_hivar_delta(1:2:end, :) = -coefEsts_out_log_hivar_delta(1:2:end, :);




%% plot

% combine
combined_prediction_strength_deltas = [subj_mevar_predstr(:);subj_hivar_predstr(:)];
combined_log_delta_coefs = [coefEsts_out_log_mevar_delta(:);coefEsts_out_log_hivar_delta(:)];
combined_normal_delta_coefs = [coefEsts_out_normal_mevar_delta(:);coefEsts_out_normal_hivar_delta(:)];

coefs_less_than_60_idx = abs(combined_log_delta_coefs)<60 & abs(combined_normal_delta_coefs)<60;

% normal figure
%
figure; hold on
[r,p] = fit_line(combined_prediction_strength_deltas(coefs_less_than_60_idx), combined_normal_delta_coefs(coefs_less_than_60_idx), 0); title norm
ylim([-60 60]); axis square

plot(xlim, [0 0], 'k--'); plot([0 0], ylim, 'k--');
xhold = subj_mevar_predstr(:); yhold = coefEsts_out_normal_mevar_delta(:);
plot(xhold(yhold<60 & yhold>-60),yhold(yhold<60 & yhold>-60),'go')
plot(xlim, nanmean(yhold(yhold<60 & yhold>-60)).*[1 1], 'g--'); plot(nanmean(xhold(yhold<60 & yhold>-60)).*[1 1], ylim, 'g--');
xhold = subj_hivar_predstr(:); yhold = coefEsts_out_normal_hivar_delta(:);
plot(xhold(yhold<60 & yhold>-60),yhold(yhold<60 & yhold>-60),'bo')
plot(xlim, nanmean(yhold(yhold<60 & yhold>-60)).*[1 1], 'b--'); plot(nanmean(xhold(yhold<60 & yhold>-60)).*[1 1], ylim, 'b--');

xlabel('Prediction strength delta')
ylabel('Coefficient delta')
title(['Normal coef; r = ' num2str(r) '; p = ' num2str(p)])
axis square

% log figure
%
figure; hold on
[r,p] = fit_line(combined_prediction_strength_deltas(coefs_less_than_60_idx), combined_log_delta_coefs(coefs_less_than_60_idx), 0); title norm
ylim([-60 60]); axis square

plot(xlim, [0 0], 'k--'); plot([0 0], ylim, 'k--');
xhold = subj_mevar_predstr(:); yhold = coefEsts_out_log_mevar_delta(:);
plot(xhold(yhold<60 & yhold>-60),yhold(yhold<60 & yhold>-60),'go')
plot(xlim, nanmean(yhold(yhold<60 & yhold>-60)).*[1 1], 'g--'); plot(nanmean(xhold(yhold<60 & yhold>-60)).*[1 1], ylim, 'g--');
xhold = subj_hivar_predstr(:); yhold = coefEsts_out_log_hivar_delta(:);
plot(xhold(yhold<60 & yhold>-60),yhold(yhold<60 & yhold>-60),'bo')
plot(xlim, nanmean(yhold(yhold<60 & yhold>-60)).*[1 1], 'b--'); plot(nanmean(xhold(yhold<60 & yhold>-60)).*[1 1], ylim, 'b--');

xlabel('Prediction strength delta')
ylabel('Coefficient delta')
title(['Normal coef; r = ' num2str(r) '; p = ' num2str(p)])
title(['Logistic coef; r = ' num2str(r) '; p = ' num2str(p)])





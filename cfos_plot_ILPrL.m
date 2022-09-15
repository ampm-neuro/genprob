
% input data
load('cfos_cts.mat')
%cfos_conditions % only need to run the first time

incl_grps = 1:6;
incl_subjs = 19:24;

%% clean density
%
%
all_density_roi = all_area_Prl;
%{
for islice = included_slices
    nnan_idx = ~isnan(all_area_acc(:,islice)) & ~isnan(all_cfos_cts_acc(:,islice));
    outlier_index = mdist_outlier([all_area_acc(nnan_idx,islice) all_cfos_cts_acc(nnan_idx,islice)], 5); title(num2str(islice));
    all_density_acc(outlier_index,islice) = nan;
    close
end
%}
%all_density_acc = all_density_summary_acc; included_slices = 1;
all_density_roi(setdiff(1:30, incl_subjs),:) = nan;



% roi s
included_slices = 1:size(all_density_roi,2);




%% parse data by group ids

log_coef_cell = cell(1,6);
norm_coef_cell = cell(1,6);
density_cell = cell(1,6);
for igrp = 1:6
    
    % density
    density_cell{igrp} = nanmean(all_density_roi(subject_condition==igrp, included_slices), 2); % drop bad slices
    %density_cell{igrp} = mean(all_density_roi(subject_condition==igrp, included_slices), 2); % drop bad subjs
    
    % log
    if ismember(igrp, [1 4])
        log_coef_cell{igrp} = -subj_coefs(subject_condition==igrp, 3); 
    else
        log_coef_cell{igrp} = subj_coefs(subject_condition==igrp, 3);
    end
    %log_coef_cell{igrp}(isnan(density_cell{igrp})) = nan;
    
    % norm
    norm_coef_cell{igrp} = subj_coefs(subject_condition==igrp, 4);
    %norm_coef_cell{igrp}(isnan(density_cell{igrp})) = nan;
    
    
end


%% plot density
figure; hold on
incl_grps_mevar = incl_grps(incl_grps<=3);
mevar_xlim = [1 4 7];
errorbar_barplot(density_cell(incl_grps_mevar), 0, mevar_xlim(ismember(1:3,incl_grps_mevar)) , 0.7.*ones(length(density_cell),3));
incl_grps_hivar = incl_grps(incl_grps>=4);
errorbar_barplot(density_cell(incl_grps_hivar), 0, mevar_xlim(ismember(4:6,incl_grps_hivar))+1, 0.7.*ones(length(density_cell),3));
xlim auto
title(num2str(included_slices))


%% plot log coefs
%{
figure; hold on; 
incl_grps_mevar = incl_grps(incl_grps<=3);
mevar_xlim = [1 4 7];
errorbar_barplot(log_coef_cell(incl_grps_mevar), 0, mevar_xlim(ismember(1:3,incl_grps_mevar)) , 0.7.*ones(length(log_coef_cell),3));
incl_grps_hivar = incl_grps(incl_grps>=4);
errorbar_barplot(log_coef_cell(incl_grps_hivar), 0, mevar_xlim(ismember(4:6,incl_grps_hivar))+1, 0.7.*ones(length(log_coef_cell),3));
xlim auto
title('log coefs')


%% plot norm coefs
figure; hold on; 
incl_grps_mevar = incl_grps(incl_grps<=3);
mevar_xlim = [1 4 7];
errorbar_barplot(norm_coef_cell(incl_grps_mevar), 0, mevar_xlim(ismember(1:3,incl_grps_mevar)) , 0.7.*ones(length(norm_coef_cell),3));
incl_grps_hivar = incl_grps(incl_grps>=4);
errorbar_barplot(norm_coef_cell(incl_grps_hivar), 0, mevar_xlim(ismember(4:6,incl_grps_hivar))+1, 0.7.*ones(length(norm_coef_cell),3));
xlim auto
title('norm coefs')
%}



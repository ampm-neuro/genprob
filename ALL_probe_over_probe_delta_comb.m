% ALL_probe_over_probe_delta_comb
% plots probe 2 probe change for both groups

[all_deltas_mevar, all_waits_mevar, all_logcell_mevar, all_normcell_mevar] = ALL_probe_over_probe_delta('train_mevar_fin');
[all_deltas_hivar, all_waits_hivar, all_logcell_hivar, all_normcell_hivar] = ALL_probe_over_probe_delta('train_hivar_fin');

figure; 

subplot(2,2, 1:2);hold on
errorbar_plot_lineonly(all_logcell_mevar(1:6), 1); 
errorbar_plot_lineonly(all_logcell_hivar(1:6), 1);
plot(xlim, [0 0], 'k--')
title log; 

subplot(2,2, 3:4); hold on
errorbar_plot_lineonly(all_normcell_mevar(1:6), 1); 
errorbar_plot_lineonly(all_normcell_hivar(1:6), 1);
plot(xlim, [0 0], 'k--')
title norm; 

%{
all_logcell_mevar_hold = all_logcell_mevar;
nan_idx = zeros(size(all_logcell_mevar_hold{1},1), size(all_logcell_mevar_hold,2)-1);
for icell=[1 6]
    nan_idx(:,icell) = isnan(all_logcell_mevar_hold{icell});
end
nan_idx = sum(nan_idx,2)>0;
for icell=[1 6]
    all_logcell_mevar_hold{icell} = all_logcell_mevar_hold{icell}(~nan_idx);
end

all_logcell_hivar_hold = all_logcell_hivar;
nan_idx = zeros(size(all_logcell_hivar_hold{1},1), size(all_logcell_hivar_hold,2)-1);
for icell=[1 6]
    nan_idx(:,icell) = isnan(all_logcell_hivar_hold{icell});
end
nan_idx = sum(nan_idx,2)>0;
for icell=[1 6]
    all_logcell_hivar_hold{icell} = all_logcell_hivar_hold{icell}(~nan_idx);
end
%}
figure; hold on
errorbar_plot(all_logcell_mevar([1 6]), 1, [], ([57 116 29]./255).*ones(size(all_logcell_mevar{1},1), 3)); 
errorbar_plot(all_logcell_hivar([1 6]), 1, [], ([41 91 173]./255).*ones(size(all_logcell_hivar{1},1), 3)); 
plot(xlim, [0 0], 'k--')
title log; 
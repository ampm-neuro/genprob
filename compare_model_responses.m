[logist_coef_cell_corrected_mevar, norm_coef_cell_mevar] = ALL_probe_vs_behavior_fit('train_mevar');
[logist_coef_cell_corrected_hivar, norm_coef_cell_hivar] = ALL_probe_vs_behavior_fit('train_hivar');

%% prep
xpos = [1 4 7 10 13];
dark_gray = .3.*[1 1 1];
light_gray = .8.*[1 1 1];
xtl = [{'Pre-train'} {'Post-problem01'} {'Post-train'} {'1 week'} {'Muted'}];




%% logistic
mevar_hold_log = logist_coef_cell_corrected_mevar([1 2 7 8 9]);
hivar_hold_log = logist_coef_cell_corrected_hivar([1 2 7 8 9]);

figure; hold on
errorbar_plot_noline(mevar_hold_log, 0, xpos, dark_gray)
errorbar_plot_noline(hivar_hold_log, 0, xpos+1, light_gray)
xlim([xpos(1)-1 xpos(end)+2])
xticks(xpos+0.5)
xticklabels(xtl)
plot(xlim, [1 1].*0, 'k--')
ylim([-15 21])
title logistic

% ttests and asterisks
for ittest = 1:length(mevar_hold_log)
    
    %ttests
    %
    [~,pval] = ttest(mevar_hold_log{ittest});
    %sig_asterisks(pval, xpos(ittest)+0, 18)
    [~,pval] = ttest(hivar_hold_log{ittest});
    %sig_asterisks(pval, xpos(ittest)+1, 18)
    %}
    
    % ttest2
    [~,pval] = ttest2(mevar_hold_log{ittest}, hivar_hold_log{ittest});
    sig_asterisks(pval, xpos(ittest)+0.25, 19)
end


%% normal
mevar_hold_norm = norm_coef_cell_mevar([1 2 7 8 9]);
hivar_hold_norm = norm_coef_cell_hivar([1 2 7 8 9]);

figure; hold on
errorbar_plot_noline(mevar_hold_norm, 0, xpos, dark_gray)
errorbar_plot_noline(hivar_hold_norm, 0, xpos+1, light_gray)
xlim([xpos(1)-1 xpos(end)+2])
xticks(xpos+0.5)
xticklabels(xtl)
plot(xlim, [1 1].*0, 'k--')
ylim([-15 23])
title normal

% ttests and asterisks
for ittest = 1:length(mevar_hold_norm)
    
    %ttests
    %
    [~,pval] = ttest(mevar_hold_norm{ittest});
    %sig_asterisks(pval, xpos(ittest)+0, 20)
    [~,pval] = ttest(hivar_hold_norm{ittest});
    %sig_asterisks(pval, xpos(ittest)+1, 20)
    %}
    
    % ttest2
    [~,pval] = ttest2(mevar_hold_norm{ittest}, hivar_hold_norm{ittest});
    sig_asterisks(pval, xpos(ittest)+0.25, 21)
    
end


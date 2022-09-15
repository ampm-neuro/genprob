
mevar_subjs = [{'651049m1'} {'658648m2'} {'691359m1'} {'690330m1'} {'691359m2'} {'696944m3'} {'699437m4'} {'699438m3'}];
hivar_subjs = [{'687034m4'} {'683472m2'} {'699437m2'} {'699438m1'} {'683472m3'}];
green_color = [0 180 0]./255;
blue_color =  [0 115 207]./255;

% load('cell_corrs_1&4&16&19_may2022.mat')
[first_sessions_mevar_early, last_sessions_mevar_early, corrs_mevar_early, subj_idx_mevar_early] = pairwise_session_remapping('mevar', mevar_subjs, [1 4]); 
[first_sessions_hivar_early, last_sessions_hivar_early, corrs_hivar_early, subj_idx_hivar_early] = pairwise_session_remapping('hivar', hivar_subjs, [1 4]);
[first_sessions_mevar_late, last_sessions_mevar_late, corrs_mevar_late, subj_idx_mevar_late] = pairwise_session_remapping('mevar', mevar_subjs, [16 19]); 
[first_sessions_hivar_late, last_sessions_hivar_late, corrs_hivar_late, subj_idx_hivar_late] = pairwise_session_remapping('hivar', hivar_subjs, [16 19]);
%[first_sessions_mevar_week, last_sessions_mevar_week, corrs_mevar_week, subj_idx_mevar_week] = pairwise_session_remapping('mevar', mevar_subjs, [19 20]); 
%[first_sessions_hivar_week, last_sessions_hivar_week, corrs_hivar_week, subj_idx_hivar_week] = pairwise_session_remapping('hivar', hivar_subjs, [19 20]);

figure; hold on
%colors = distinguishable_colors(length(subject_ids));
%for isubj = 1:length(subject_ids)
%    plot(atanh(cell_corrs(subj_idx==isubj)), 'o', 'color', colors(isubj,:));
%    plot(mean(atanh(cell_corrs(subj_idx==isubj))), '.', 'color', colors(isubj,:), 'markersize', 40)
%end
legend_trick([green_color; blue_color], '-')
errorbar_plot_noline([{atanh(corrs_hivar_early)} {atanh(corrs_hivar_late)}], 0, [1 2], blue_color+.1, blue_color+.1)
errorbar_plot_noline([{atanh(corrs_mevar_early)} {atanh(corrs_mevar_late)}], 0, [1 2], green_color+.1, green_color+.1)
errorbar_plot_lineonly([{atanh(corrs_hivar_early)} {atanh(corrs_hivar_late)}], 0, [1 2], blue_color, blue_color)
errorbar_plot_lineonly([{atanh(corrs_mevar_early)} {atanh(corrs_mevar_late)}], 0, [1 2], green_color, green_color)
tic_vect = [-.99 -.95 -.9 -.8 -.6 -.3 0 .3 .6 .8 .9 .95 .99 .999]; 
ylim(atanh([-.6 .9])); yticks(atanh(tic_vect)); yticklabels(tic_vect); ylim(atanh([-.95 .995]))
xticks([1 2]); xlim([.5 2.5])
hold on; plot(xlim, [0 0], 'k--')
xticklabels([{'Problem 1'}, {'Problem 6'}])
xlabel('Before and after probe comparison')
ylabel('Similarity of neuronal activity (r)')
title('Stability of neuronal activity after learning')
legend({'Predictable', 'Unpredictable'})
legend('location', 'northeastoutside')
function corr_train_probe(probe_num, day_delay, new_learn_num)
% correlates number of training days with probe dvs

% probe dvs
[mus, sigmas, multipliers, intercepts, gofs, subjects_probe] = plot_probes_means(probe_num, day_delay, new_learn_num, 0);

% number of learning days
[nld, subjects_nld] = num_learning_days;
nld = nld(ismember(subjects_nld, subjects_probe),:);
nld_sum = sum(nld,2);

% plot
figure;
subplot(3,2,1); [R, p] = fit_line(nld_sum, mus); xlabel('Total training days')
title(['r=' num2str(R) ', p=' num2str(p)]); ylabel('Probe mus')
set(gca,'TickLength',[0, 0]); box off; axis square;
subplot(3,2,2); [R, p] = fit_line(nld_sum, sigmas); xlabel('Total training days')
title(['r=' num2str(R) ', p=' num2str(p)]); ylabel('Probe sigmas')
set(gca,'TickLength',[0, 0]); box off; axis square;
subplot(3,2,3); [R, p] = fit_line(nld_sum, multipliers); xlabel('Total training days')
title(['r=' num2str(R) ', p=' num2str(p)]); ylabel('Probe multipliers')
set(gca,'TickLength',[0, 0]); box off; axis square;
subplot(3,2,4); [R, p] = fit_line(nld_sum, intercepts); xlabel('Total training days')
title(['r=' num2str(R) ', p=' num2str(p)]); ylabel('Probe intercepts')
set(gca,'TickLength',[0, 0]); box off; axis square;
subplot(3,2,5); [R, p] = fit_line(nld_sum, gofs); xlabel('Total training days')
title(['r=' num2str(R) ', p=' num2str(p)]); ylabel('Probe gofs')
set(gca,'TickLength',[0, 0]); box off; axis square;
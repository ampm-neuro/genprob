[wait_durations, wd_freq, freq_numbers] = wait_times_prep(trl_mtx,2); 
[fit_yvals, coefEsts, modelFun] = ampm_normal_logistic_fit_algo(freq_numbers, wait_durations);
figure; hold on
plot(freq_numbers, wait_durations,'o'); 
plot((1:0.001:41), modelFun(coefEsts, (1:0.001:41)'));
set(gca,'TickLength',[0, 0]); box off;
coefEsts
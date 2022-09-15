
% plot trials
plot_trials_trlmtx(trl_mtx) %probe
plot_trials_trlmtx_altcolor(trl_mtx) % problem

% plot single session
figure; wait_times_plot_opto(trl_mtx,1);

% plot all probes smooth, and coefficient changes
all_probe_smooth_meta


% plot learning curves during a single super session
ALL_wait_times_sesh_plot_norm({'opto-01'}); % can select sessions near top of code

% problem performance
 ALL_opto_discrim_folder([1 6], 'train_hivar_optoExp_hpc'); % Input desired session numbers

% probe over probe change
[probe_brON, probe_bpON, probe_brOFF, probe_bpOFF, probe_subj_ids, probe_problemNums] = probe_change_opto(1:7);

% compare probe behavior and problem behavior with opto (eg if laser affects
% last day problem performance, does it affect probe performance?)
probe_vs_behavior_opto


% plot average probe responses with and without laser ON
current_folder = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_all_optoExp_hpc\';
for iprobe = [1 2]
    figure; 
    if iprobe<7
        [all_fixed_effects, all_beta_pvals] = ...
        plot_allprobes_opto(get_file_paths_targeted(current_folder, {'preprobe_opto', ['0' num2str(iprobe) '_']}));
    elseif iprobe<9
        [all_fixed_effects, all_beta_pvals] = ...
        plot_allprobes_opto(get_file_paths_targeted(current_folder, {'postprobe_opto', ['0' num2str(iprobe-6) '_']}));
    elseif iprobe==9
        %[all_fixed_effects, all_beta_pvals] = ...
        %plot_allprobes_opto(get_file_paths_targeted(current_folder, {'postprobe_opto', ['0' num2str(iprobe-6) '_']}));
    end
title(['Probe 0' num2str(iprobe)])    
end

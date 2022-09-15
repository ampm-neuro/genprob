function hm_cells = tw_activity_plot_trial_hm_full_warp(trl_mtx, trl_idx, frame_times, traces, neurons, trials, tw, tses)
% a subplot of all tw_activity_plot_trial_hm alignments for each cell
% input tses is time series event spacing used to warp traces

hm_cells = cell(1, length(neurons));

for ineuron = 1:length(neurons)
    
    figure; 
    hm_cells{ineuron} = tw_activity_plot_trial_hm_warp(trl_mtx, trl_idx, frame_times, traces, neurons(ineuron), trials, tw, tses);
    caxis([0 1])
end
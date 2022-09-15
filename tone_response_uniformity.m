
% quantify similarity of responses across tone frequencies for each neuron
tses = [0.2 1.1 1.0 2.0 8.5 2.0];
all_matrices = tw_activity_plot_tone_hm_warp(trl_mtx, trl_idx, medass_cell, frame_times, traces,1:size(traces,1), 1:size(trl_mtx,1), 4, [5 25], tses);
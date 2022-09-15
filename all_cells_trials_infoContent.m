function all_info_content = all_cells_trials_infoContent(neurons, session_mtx_cell)
% compute and plot average infocontent

% first 41 trials from each session

% load('image_pop_mds_full_save.mat', 'session_mtx_cell')

all_info_content = [];

for icell = neurons

    all_info_content = [all_info_content; ...
        all_trials_one_cell_infoContent(icell, session_mtx_cell)'];
    
end


figure; hold on
plot(nanmean(all_info_content), 'linewidth', 2, 'color', [0 0 0])
se = nanstd(all_info_content)./sqrt(sum(~isnan(all_info_content)));
plot(nanmean(all_info_content) + se, 'linewidth', 1, 'color', 0.7.*[1 1 1])
plot(nanmean(all_info_content) - se, 'linewidth', 1, 'color', 0.7.*[1 1 1])


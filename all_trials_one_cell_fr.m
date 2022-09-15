function [fullcell, fullcell_norm] = all_trials_one_cell_fr(neurons, session_mtx_cell, tses)
% load session_mtx_cell computed inside image_pop_mds 

% timewarp info
tses = [0 tses];
tses = tses.*100; % 100 frames per second

% number of trials in each session
session_trials = nan(length(session_mtx_cell),1);
for isesh = 1:length(session_mtx_cell)
    session_trials(isesh) = size(session_mtx_cell{isesh}{1},1);
end
cum_session_trials = cumsum(session_trials);

% iterate through neurons
for cell_num = neurons
    fullcell = []; 
    
    for isesh = 1:length(session_mtx_cell)
        %fullcell = [fullcell; session_mtx_cell{isesh}{cell_num}(1:41,:)];
        fullcell = [fullcell; session_mtx_cell{isesh}{cell_num}];
    end
    
    figure
    hold on
    
    % plot rates
    fullcell_norm = norm_mtx(fullcell')';
    
    imagesc(fullcell_norm)
    %ampm_pcolor(fullcell_norm)
    set(gca, 'YDir','reverse')
    %ampm_pcolor(fullcell)
    axis([1 size(fullcell_norm,2) 1 size(fullcell_norm,1)])
    
    % resize
    %set(gcf, 'Position', [1000 380 560 958])
    
    % plot trial bounds
    %
    for isesh = 2:length(session_trials)
        plot(xlim, cum_session_trials(isesh-1).*[1 1], 'k-', 'linewidth', 1.1)
    end
    
    % plot red lines
    ctses = cumsum(tses);
    for rl = ctses([1 2 3 4 6])
        plot(rl.*[1 1], ylim, 'r-', 'linewidth', 2)
    end
    plot(ctses(5).*[1 1], ylim, 'r--', 'linewidth', 2)
    
    % axis labels
    xticks(ctses)
    xticklabels_universe = {'NP1', 'NP2', 'HE1', 'TO', 'RD', 'RWD', 'TE'};
    xticklabels(xticklabels_universe)
    yticks(mean([cum_session_trials'; [0 cum_session_trials(1:end-1)']]))
    yticklabels(1:length(session_trials))
    ylabel('Session')
    
    
end
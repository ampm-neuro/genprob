function tw_activity_plot_trial_hm_full(trl_mtx, trl_idx, medass_cell, frame_times, traces, neurons, trials, tw)
% a subplot of all tw_activity_plot_trial_hm alignments for each cell

% str
str_input = {'NP onset','NP offset','Tone on','Head entry',...
    'Reward delivery','Reward receipt','Trial end'};

for ineuron = 1:length(neurons)
    
    figure; 
    
    alignments = [1 2 3 4 5 7];
    subplot_order = [1 2 4 3 5 6];
    for ialign = 1:length(alignments)
    
        subplot(1,length(alignments), subplot_order(ialign))
        hold on
        tw_activity_plot_trial_hm(trl_mtx, trl_idx, medass_cell, frame_times, traces, neurons(ineuron), trials, alignments(ialign), tw)
        
        if subplot_order(ialign) ~= 1
            yticks([])
            ylabel([])
        end
        
        title(str_input(alignments(ialign)))
        set(gcf,'Position', [78        1037        2412         301])
        caxis([0 1])
    end
end
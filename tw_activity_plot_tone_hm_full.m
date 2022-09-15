function tw_activity_plot_tone_hm_full(trl_mtx, medass_cell, frame_times, traces, neurons, trials, tw)
% a subplot of all tw_activity_plot_tone_hm alignments for each cell

% str
str_input = {'NP onset','NP offset','Tone on','Head entry',...
    'Reward delivery','Reward receipt','Trial end'};


for ineuron = 1:length(neurons)
    
    figure; 
    
    for ialign = 1:7
    
        subplot(1,7, ialign)
        hold on
        tw_activity_plot_tone_hm(trl_mtx, medass_cell, frame_times, traces, neurons(ineuron), trials, ialign, tw)
        
        yticks([])
        %xticks([])
        ylabel([])
        %xlabel([])
        colorbar off
        
        title(str_input(ialign))
        set(gcf,'Position', [78        1037        2412         301])
    end
end
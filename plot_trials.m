function plot_trials(medass_cell)
% plots major events of each trial as a row

% rich tones
load('unqfrq41.mat', 'unqfrq41');
rich_tones = unqfrq41(16:22);
clearvars unqfrq41;

figure; hold on; set(gca,'TickLength',[0, 0]);
colors = get(gca,'ColorOrder');
max_trial_length = [];
valid_trial_count = 0;

%completed trials only
medass_cell{2} = medass_cell{2}(1:length(medass_cell{4}));

%iterate through trials
for itrl = 1:length(medass_cell{2})
    %trial time bounds
    trial_start_np = max(medass_cell{7}(medass_cell{7} <= medass_cell{6}(itrl))); %nose poke
    trial_start_he = max(medass_cell{12}(medass_cell{12} > trial_start_np & medass_cell{12} < medass_cell{4}(itrl))); %head entry
    %trial_start_he = min(medass_cell{12}(medass_cell{12} > trial_start_np)); %head entry
    trial_start_fixed = trial_start_he + medass_cell{1}(6)/100;
    
    % for sessions where tone came on after head entry
    %{
    if length(medass_cell{1})>=32 && medass_cell{1}(32)>0
       
        trial_tone_jitter = medass_cell{9}(itrl) - trial_start_he;
        tsf = trial_start_fixed;
        trial_start_fixed = trial_start_fixed + medass_cell{1}(32) + trial_tone_jitter;
        
        hold_temp = [hold_temp trial_start_fixed - tsf]
        
    end
    %}
    
    trial_end = medass_cell{4}(itrl);
    trial_end_he_offset = medass_cell{13}(medass_cell{13}>trial_start_fixed & medass_cell{13}<medass_cell{4}(itrl))';
    trial_end = min([trial_end trial_end_he_offset]);    
    
    if trial_end < trial_start_fixed
        continue
    end
    
    %xlim details
    max_trial_length = max([max_trial_length (trial_end-trial_start_fixed)]);
    
    %head entrys
    head_entry = medass_cell{12}(trl_time(medass_cell{12}, trial_start_he, trial_end))';
    
    if ~isempty(head_entry)
        valid_trial_count = valid_trial_count+1;
        
        %head entry and exits
        head_exit = medass_cell{13}(trl_time(medass_cell{13}, trial_start_he, trial_end))';
        head_entry = head_entry - trial_start_fixed;
        plot([head_entry;head_entry], repmat([.5;1.5], 1, length(head_entry))+(valid_trial_count-1), 'Color', colors(2,:))
        if ~isempty(head_exit)
            head_exit = head_exit - trial_start_fixed;
        else
            head_exit = trial_end - trial_start_fixed;
        end
        plot([head_exit;head_exit], repmat([.5;1.5], 1, length(head_exit))+(valid_trial_count-1), 'Color', [0 0 0], 'linewidth', 1)
        
        %head entry line
        trial_head_entry = min(head_entry);
        trial_head_exit = min(head_exit(head_exit>trial_head_entry));
        
        %nose pokes
        nose_poke_onset_times = medass_cell{7}(trl_time(medass_cell{7}, trial_start_np, trial_end))';
        nose_poke_onset_times = nose_poke_onset_times - trial_start_fixed;
        nose_poke_offset_times = medass_cell{8}(trl_time(medass_cell{8}, trial_start_np, trial_end))';
        nose_poke_offset_times = nose_poke_offset_times - trial_start_fixed;
        plot([nose_poke_onset_times;nose_poke_onset_times], repmat([0.5;1.5], 1, length(nose_poke_onset_times))+(valid_trial_count-1), 'Color', colors(1,:))
        plot([nose_poke_offset_times;nose_poke_offset_times], repmat([0.5;1.5], 1, length(nose_poke_onset_times))+(valid_trial_count-1), 'Color', colors(1,:))
        plot([nose_poke_onset_times;nose_poke_offset_times], repmat([1.0;1.0], 1, length(nose_poke_onset_times))+(valid_trial_count-1), 'Color', colors(1,:))

        %rewards and waits
        rewards = medass_cell{14}(trl_time(medass_cell{14}, trial_start_he, trial_end))';
        rewards = rewards - trial_start_fixed;
        if ~isempty(rewards)
            %plot([trial_head_entry; trial_head_exit], [1; 1]+(valid_trial_count-1), 'Color', 0.7.*[1 1 1])
            if ismember(floor(medass_cell{10}(itrl)), rich_tones)
                plot([trial_head_entry; trial_head_exit], [1; 1]+(valid_trial_count-1), 'Color', [174 116 116]./255, 'linewidth', 1)
            else
                plot([trial_head_entry; trial_head_exit], [1; 1]+(valid_trial_count-1), 'Color', 0.5.*[1 1 1], 'linewidth', 1)
            end
            plot([rewards;rewards], repmat([0.5;1.5], 1, length(rewards))+(valid_trial_count-1), 'Color', [22 117 15]./255, 'linewidth', 2.5)
        elseif medass_cell{3}(itrl) == 0
            if ismember(floor(medass_cell{10}(itrl)), rich_tones)
                plot([trial_head_entry; trial_head_exit], [1; 1]+(valid_trial_count-1), 'Color', [166 28 28]./255, 'linewidth', 1.5)
            else
                plot([trial_head_entry; trial_head_exit], [1; 1]+(valid_trial_count-1), 'Color', 0.0.*[1 1 1], 'linewidth', 1.5)
            end
        else
            %plot([trial_head_entry; trial_head_exit], [1; 1]+(valid_trial_count-1), 'Color', 0.5.*[1 1 1])
            if ismember(floor(medass_cell{10}(itrl)), rich_tones)
                plot([trial_head_entry; trial_head_exit], [1; 1]+(valid_trial_count-1), 'Color', [174 116 116]./255, 'linewidth', 1)
            else
                plot([trial_head_entry; trial_head_exit], [1; 1]+(valid_trial_count-1), 'Color', 0.5.*[1 1 1], 'linewidth', 1)
            end
        end
        
        
        % opto patch
        if sum(medass_cell{17}>0) && sum(medass_cell{17})<(length(medass_cell{17})/2) % if opto session
            if medass_cell{17}(itrl)==1 % if opto trial
                patch([min(nose_poke_onset_times) trial_head_exit trial_head_exit min(nose_poke_onset_times)],[0.5 0.5 1.5 1.5]+(valid_trial_count-1), [47 141 255]./255, 'FaceAlpha', .3)
            end
        end
        
        
    end

end

%aesthetics
ylim([0.5 valid_trial_count+.5])
hold on; plot([0;0]-medass_cell{1}(6), ylim, '-', 'color', [.7 .7 .7]) %head entry
hold on; plot([0;0], ylim, '-', 'color', [.7 .7 .7]) %fixed delay    
xlim([-6 max_trial_length+1])
xlabel('Time (s)')
ylabel('Trial')


end




function trl_time_out = trl_time(cell_of_interest, trl_start, trl_end)
    trl_time_out = cell_of_interest>=trl_start & cell_of_interest<=trl_end;
end
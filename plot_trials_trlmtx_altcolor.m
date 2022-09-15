function plot_trials_trlmtx_altcolor(trl_mtx)
% plots major events of each trial as a row

% rich tones
load('unqfrq41.mat', 'unqfrq41');
rich_tones = unqfrq41(13:29);
clearvars unqfrq41;

% figure prep
figure; hold on; set(gca,'TickLength',[0, 0]);
colors = get(gca,'ColorOrder');
valid_trial_count = 0;
max_trial_length = [];


%iterate through trials
for itrl = 1:size(trl_mtx,1)
    valid_trial_count = valid_trial_count+1;
    
    %trial time bounds
    trial_start_np = sum(trl_mtx(itrl, [1 6])); % nose poke ON
    np_offset = sum(trl_mtx(itrl, [1 9])); % nose poke OFF
    head_entry = sum(trl_mtx(itrl, [1 10])); %head entry ON
    tone_on = sum(trl_mtx(itrl, [1 7])); % tone on
    rand_del_start = trl_mtx(itrl, 1);
    rwd_delivery = sum(trl_mtx(itrl, [1 11])); % reward delivery
    trl_end = sum(trl_mtx(itrl, [1 12])); % head entry OFF or post rwd
    
    % align to head_entry
    trial_start_np = trial_start_np - head_entry;
    np_offset = np_offset - head_entry;
    tone_on = tone_on - head_entry;
    rwd_delivery = rwd_delivery - head_entry;
    trl_end = trl_end - head_entry;
    rand_del_start = rand_del_start - head_entry;
    head_entry = head_entry - head_entry;

    % plot events 
    plot([head_entry;head_entry], repmat([.5;1.5], 1, length(head_entry))+(valid_trial_count-1), 'Color', colors(2,:))
    plot([trl_end;trl_end], repmat([.5;1.5], 1, length(trl_end))+(valid_trial_count-1), 'Color', [0 0 0], 'linewidth', 1)
    plot([trial_start_np;trial_start_np], repmat([0.5;1.5], 1, length(trial_start_np))+(valid_trial_count-1), 'Color', colors(1,:))
    plot([trial_start_np;trial_start_np], repmat([1.0;1.0], 1, length(trial_start_np))+(valid_trial_count-1), 'Color', colors(1,:))
    plot([np_offset;np_offset], repmat([0.5;1.5], 1, length(np_offset))+(valid_trial_count-1), 'Color', colors(1,:))
    plot([trial_start_np;np_offset], repmat([1.0;1.0], 1, length(trial_start_np))+(valid_trial_count-1), 'Color', colors(1,:))
    plot([tone_on;tone_on], repmat([.5;1.5], 1, length(tone_on))+(valid_trial_count-1), 'Color', colors(2,:))
    plot([rand_del_start;rand_del_start], repmat([.5;1.5], 1, length(rand_del_start))+(valid_trial_count-1), 'g')
    
        
    %rewards and waits colored by trial type
    if ~isnan(trl_mtx(itrl,11)) % rewarded trial
        if ismember(floor(trl_mtx(itrl,2)), rich_tones)
            plot([head_entry; trl_end], [1; 1]+(valid_trial_count-1), 'Color', [141 2 31]./255, 'linewidth', 1)
        else
            plot([head_entry; trl_end], [1; 1]+(valid_trial_count-1), 'Color', 0.5.*[1 1 1], 'linewidth', 1)
        end
        plot([rwd_delivery;rwd_delivery], repmat([0.5;1.5], 1, length(rwd_delivery))+(valid_trial_count-1), 'Color', [22 117 15]./255, 'linewidth', 2.5)
    elseif trl_mtx(itrl,3) == 0 % reward unavailable
        if ismember(floor(trl_mtx(itrl,2)), rich_tones)
            plot([head_entry; trl_end], [1; 1]+(valid_trial_count-1), 'Color', [141 2 31]./255, 'linewidth', 2)
        else
            plot([head_entry; trl_end], [1; 1]+(valid_trial_count-1), 'Color', 0.0.*[1 1 1], 'linewidth', 2)
        end
    else % abandoned trial
        if ismember(floor(trl_mtx(itrl,2)), rich_tones)
            plot([head_entry; trl_end], [1; 1]+(valid_trial_count-1), 'Color', [141 2 31]./255, 'linewidth', 1)
        else
            plot([head_entry; trl_end], [1; 1]+(valid_trial_count-1), 'Color', 0.5.*[1 1 1], 'linewidth', 1)
        end
    end
        
        
    % opto patch
    %
    if sum(trl_mtx(:,13))>0 && sum(trl_mtx(:,13))<(size(trl_mtx,1)/2) % if opto session
        if trl_mtx(itrl,13)==1 % if opto trial
            patch([tone_on-0.2 trl_end trl_end tone_on-0.2],[0.5 0.5 1.5 1.5]+(valid_trial_count-1), [47 141 255]./255, 'FaceAlpha', .3)
        end
    end
    %}
        
    
    %xlim details
    max_trial_length = max([max_trial_length (trl_end-trial_start_np)]);
        
    end



%aesthetics
ylim([0.5 valid_trial_count+.5])
%hold on; plot([0;0]-medass_cell{1}(6), ylim, '-', 'color', [.7 .7 .7]) %head entry
hold on; plot([0;0], ylim, '-', 'color', [.7 .7 .7]) %fixed delay    
xlim([-6 max_trial_length+1])
xlabel('Time (s)')
ylabel('Trial')


end




function trl_time_out = trl_time(cell_of_interest, trl_start, trl_end)
    trl_time_out = cell_of_interest>=trl_start & cell_of_interest<=trl_end;
end
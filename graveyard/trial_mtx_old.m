function trl_mtx = trial_mtx(medass_cell)
%creates a matrix (trials, details)
%
% OUTPUT COLUMNS
%
% Trial details:
% col 1 = trial start time (end of fixed delay; raw time)
% col 2 = tone frequency and volume (freq.vol)
% col 3 = reward available (0 = no; 1 = yes)
% col 4 = delay duration
% col 5 = licks during fixed delay
% col 13 = opto (0 = OFF; 1 = ON)
%
% Trial event times (relative to trial start):
% col 6 = nose poke onset
% col 7 = tone on
% col 8 = tone off
% col 9 = nose poke offset
% col 10 = head entry
% col 11 = reward (or nan)
% col 12 = trial end (or withdrawal)

%valid and completed trials only
try
    medass_cell{2} = medass_cell{2}(1:length(medass_cell{4}));
catch
    %medass_cell{4};
end

%soft preallocate trial matrix
trl_mtx = nan(length(medass_cell{2}),12);

%iterate through trials
for itrl = 1:length(medass_cell{2})
    
    %trial time bounds
    trial_start_np = max(medass_cell{7}(medass_cell{7} <= medass_cell{6}(itrl))); %nose poke
    trial_start_he = min(medass_cell{12}(medass_cell{12} > trial_start_np)); %head entry
    trial_start_fixed = trial_start_he + medass_cell{1}(6)/100;
    trial_end = medass_cell{4}(itrl);
    trial_end_he_offset = medass_cell{13}(medass_cell{13}>trial_start_fixed & medass_cell{13}<medass_cell{4}(itrl))';
    trial_end = min([trial_end trial_end_he_offset]);    
    
    if trial_end < trial_start_fixed
        continue
    end
    
    %head entrys
    head_entry = medass_cell{12}(trl_time(medass_cell{12}, trial_start_he, trial_end))';
    
    if ~isempty(head_entry)
        
        %head entry and exits
        head_exit = medass_cell{13}(trl_time(medass_cell{13}, trial_start_he, trial_end))';
        head_entry = head_entry - trial_start_fixed;
        head_exit = head_exit - trial_start_fixed;
        
        %head entry line
        trial_head_entry = min(head_entry);
        trial_head_exit = min(head_exit(head_exit>trial_head_entry));
        
        %nose pokes
        nose_poke_onset_times = medass_cell{7}(trl_time(medass_cell{7}, trial_start_np, trial_end))';
        nose_poke_onset_times = nose_poke_onset_times - trial_start_fixed;
        nose_poke_offset_times = medass_cell{8}(trl_time(medass_cell{8}, trial_start_np, trial_end))';
        nose_poke_offset_times = nose_poke_offset_times - trial_start_fixed;

        %tone onset
        tone_onset_times = medass_cell{9}(trl_time(medass_cell{9}, trial_start_np, trial_end))';
        tone_onset_times = tone_onset_times - trial_start_fixed;
        tone_offset_times = medass_cell{16}(trl_time(medass_cell{16}, trial_start_np, trial_end))';
        tone_offset_times = tone_offset_times - trial_start_fixed;
        if isempty(tone_offset_times)
            tone_offset_times = trial_end - trial_start_fixed;
        end
        if length(tone_offset_times>=2)
            tone_offset_times = max(tone_offset_times);
        end
        %tone_onset_times = tone_offset_times - 0.25;
        %licks
        licks = medass_cell{15}(trl_time(medass_cell{15}, trial_start_he, trial_end))';
        licks = licks - trial_start_fixed;

        %rewards
        rewards = medass_cell{14}(trl_time(medass_cell{14}, trial_start_he, trial_end))';
        rewards = rewards - trial_start_fixed;
        rwd_dur = medass_cell{1}(3)/100;
    end
    
    
    
    %load trial matrix
    trl_mtx(itrl, 1) = trial_start_fixed; %load start (raw)
    if medass_cell{10}~=0
        trl_mtx(itrl, 2) = medass_cell{10}(itrl); %load tone
    end
    trl_mtx(itrl, 3) = medass_cell{3}(itrl); %reward availability
    trl_mtx(itrl, 4) = medass_cell{11}(itrl); %delay duration
    trl_mtx(itrl, 5) = length(medass_cell{15}(trl_time(medass_cell{15},...
            trial_start_he, trial_start_fixed)))./(trial_start_fixed-trial_start_he); %licks during fixed delay
    trl_mtx(itrl, 6) = min(nose_poke_onset_times); %NP start
    if ~isempty(tone_onset_times) && ~isempty(tone_offset_times)
        trl_mtx(itrl, 7) = tone_onset_times; %tone on
        trl_mtx(itrl, 8) = tone_offset_times; %tone end
    end
    trl_mtx(itrl, 9) = min(nose_poke_offset_times); %NP end
    trl_mtx(itrl, 10) = trial_head_entry; %head entry
    if ~isempty(rewards)
        trl_mtx(itrl, 11) = rewards; %reward delivery time
    else
        trl_mtx(itrl, 11) = nan; %no reward
    end
    trl_mtx(itrl, 12) = trial_end - trial_start_fixed; %trial end
    if max(size(medass_cell{17})) > 1
        trl_mtx(itrl, 13) = medass_cell{17}(itrl); %TTL
    end

        
end

end


function trl_time_out = trl_time(cell_of_interest, trl_start, trl_end)
    trl_time_out = cell_of_interest>=trl_start & cell_of_interest<=trl_end;
end







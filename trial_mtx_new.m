function [trl_mtx, warning_logical] = trial_mtx_new(medass_cell, varargin)
%creates a matrix (trials, details)
%
% OUTPUT COLUMNS
%
% Trial details:
% col 1 = trial start time (end of fixed delay; raw time)
% col 2 = tone frequency and volume (freq.vol)
% col 3 = reward available (0 = no; 1 = yes)
% col 4 = delay duration
% col 5 = lick rate during fixed delay (hz)
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

% path
if nargin > 1
    ppath = varargin{1};
else 
    ppath = 'string place holder';
end

% legacy path exception
lpe = 'E:\Projects\InProgress\GenProb\data\two_tone\train_mevar_imaging\651049m1\LED_gen11_mevar05\02\!2020-02-06_14h45m.Subject 651049m1';

%valid and completed trials only
try
    medass_cell{2} = medass_cell{2}(1:length(medass_cell{4}));
catch
end

%soft preallocate
trl_mtx = nan(length(medass_cell{2}),12);
fixed_delay_duration = medass_cell{1}(6)/100;

%iterate through trials
for itrl = 1:length(medass_cell{2})

        %trial time bounds
        trial_start_np = medass_cell{2}(itrl); %nose poke
        trial_end = medass_cell{4}(itrl);
        
        %if using head entry + tone onset jitter
        if length(medass_cell{1})>=32 && medass_cell{1}(32)>0 && ~contains(ppath, lpe)
            
            try
                % trial start head-entry
                trial_start_he = medass_cell{6}(itrl);

                %trial tone jitter
                trial_tone_jitter = medass_cell{9}(itrl) - trial_start_he + fixed_delay_duration;

                %start of random delay 
                trial_start_fixed = medass_cell{5}(itrl);
            
            catch
                 % trial start head entry is the first he after trial start
                trial_start_he = min(medass_cell{12}(medass_cell{12} > trial_start_np)); %head entry

                %start of random delay 
                last_he = max(medass_cell{12}(trl_time(medass_cell{12}, trial_start_np,...
                        trial_end - fixed_delay_duration)));%last head entry
                trial_start_fixed = last_he + fixed_delay_duration;
                
            end
            
        else
            
            % trial start head entry is the first he after trial start
            trial_start_he = min(medass_cell{12}(medass_cell{12} > trial_start_np)); %head entry

            %start of random delay 
            last_he = max(medass_cell{12}(trl_time(medass_cell{12}, trial_start_np,...
                        trial_end - fixed_delay_duration)));%last head entry
            trial_start_fixed = last_he + fixed_delay_duration;
        end
    
    if trial_end < trial_start_fixed
        continue
    end
    
    % reponse times
    head_entry = medass_cell{12}(trl_time(medass_cell{12}, trial_start_he, trial_end))';
        if isempty(head_entry) 
            rtime = medass_cell{14}(trl_time(medass_cell{14}, trial_start_he, trial_end));
            head_entry = rtime - medass_cell{20}(itrl) - fixed_delay_duration;
        end
        if isempty(head_entry) 
            head_entry = trial_end - medass_cell{20}(itrl) - fixed_delay_duration - 2;
        end
    if ~isempty(head_entry)
        
        % head entry
        trial_head_entry = trial_start_he - trial_start_fixed; %min(head_entry);
        
        % nose pokes
        nose_poke_onset_times = medass_cell{7}(trl_time(medass_cell{7}, trial_start_np, trial_end))';
        if isempty(nose_poke_onset_times);nose_poke_onset_times = trial_start_np; end
        nose_poke_onset_times = nose_poke_onset_times - trial_start_fixed;
        nose_poke_offset_times = medass_cell{8}(trl_time(medass_cell{8}, trial_start_np, trial_end))';
        nose_poke_offset_times = nose_poke_offset_times - trial_start_fixed;

        % tone onset
        tone_onset_times = medass_cell{9}(trl_time(medass_cell{9}, trial_start_np, trial_end))';
        tone_onset_times = tone_onset_times - trial_start_fixed;
        tone_offset_times = medass_cell{16}(trl_time(medass_cell{16}, trial_start_np, trial_end))';
        tone_offset_times = tone_offset_times - trial_start_fixed;
        if isempty(tone_offset_times)
            tone_offset_times = trial_end - trial_start_fixed;
        end
        if length(tone_offset_times)>=2
            tone_offset_times = max(tone_offset_times);
        end

        % reward times
        rewards = medass_cell{14}(trl_time(medass_cell{14}, trial_start_he, trial_end))';
        rewards = rewards - trial_start_fixed;
        if isempty(rewards) & medass_cell{3}(itrl)==1 %double check
            wait_time = (trial_end - trial_start_fixed);
            rwd_time = sum([medass_cell{11}(itrl) trial_head_entry]);
            if wait_time > rwd_time % waited long enough for reward
                rewards = rwd_time; % assign reward time
            end
        end

    end
    
    % load trial matrix    
    trl_mtx(itrl, 1) = trial_start_fixed; %load start (raw)
    if medass_cell{10}~=0
        trl_mtx(itrl, 2) = medass_cell{10}(itrl); %load tone
    end
    trl_mtx(itrl, 3) = medass_cell{3}(itrl); %reward availability
    trl_mtx(itrl, 4) = medass_cell{11}(itrl); %delay duration
    trl_mtx(itrl, 5) = length(medass_cell{15}(trl_time(medass_cell{15},...
            trial_start_he, trial_start_fixed)))./(trial_start_fixed-trial_start_he); %lick rate during fixed delay
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
        
        try
            trl_mtx(itrl, 13) = medass_cell{17}(itrl); %TTL (ONLY INTERPRETABLE FOR OPTO SESSIONS)
        catch
        end
    end
        
end

%% edit TTL if undifferentiated
warning_logical = 0;
if size(trl_mtx,2)>=13
    if any(trl_mtx(:,13)>0) && all(ismember(trl_mtx(:,13), trl_mtx(1,13))) && ~contains(ppath,'LED')
        opto_tones=[5249,5786,6377,7028,7747,8538,9411,10372,11432,12601,13888,15307,16872,18596,20496,22590,24899,27443,30247,33338];
        opto_tone_idx = ismember(floor(trl_mtx(:,2)), opto_tones);
        trl_mtx(opto_tone_idx,13) = 1;
        trl_mtx(~opto_tone_idx,13) = 0;
        disp('tones used to determine opto trials')
        warning_logical = 1;
    end
end

end




function trl_time_out = trl_time(cell_of_interest, trl_start, trl_end)
    trl_time_out = cell_of_interest>=trl_start & cell_of_interest<=trl_end;
end







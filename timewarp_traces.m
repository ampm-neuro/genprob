function [warp_traces, warp_frame_times, warp_trl_idx, warp_trl_mtx] = timewarp_traces(traces, frame_times, trl_idx, trl_mtx, time_series_event_spacing)
% rework inputs to conform to a single trial time series defined by
% 'time_series' input
%
% time time_series_event_spacing input: [nose_poke offset, head_entry, tone
% on, reward, trial end]

% catch legacy tone-on event coding
if trl_mtx(1,7)<trl_mtx(1,10)
    display('Tone-on time corrected')
    trl_mtx = correct_notone(trl_mtx);
end

% set trl_mtx(:,1) to zeros so that all timing is local
trl_mtx_1_hold = trl_mtx(:,1);
trl_mtx(:,1) = 0;

% prototypical trial event timing
trial_length_time = sum(time_series_event_spacing);
num_traces_per_trial = length(0:0.01:trial_length_time);

% preallocate output
warp_trl_mtx = trl_mtx;
warp_traces = []; 
warp_frame_times = []; 
warp_trl_idx = []; 

% trl mtx columns corresponding to time_series_event_spacing inputs
trl_mtx_col_seq = [6 9 10 7 1 11 12];

% number of frames in each new chunk according to time_series_event_spacing
num_new_chunk_frames = time_series_event_spacing./0.01;

% iterate through each trial
for itrl = 1:size(trl_mtx,1)
    if sum(trl_idx==itrl)==0
       continue 
    end
    
    % preallocate trial frames and traces
    warp_frame_times_trl = [];
    warp_traces_trl = [];
        
    % preprocess trial frames and traces
    trialft = frame_times(trl_idx == itrl);
    trial_traces = traces(:, trl_idx == itrl);
    
    % iterate through event chunks
    for itsq = 2:length(trl_mtx_col_seq)
        
        % num_new_chunk_frames
        nncf = num_new_chunk_frames(itsq-1);
        
        % update trial matrix event times (no reward)
        if trl_mtx_col_seq(itsq)==11 && isnan(warp_trl_mtx(itrl, 11))

            % rwd is nan, set trial end
            % how does trial end for a a nonrwd trial (HE withdrawl)
            % compare to trial end for a rwd trial (first lick at least 
            % 2.1s after reward pump onset (rwd pump stays on 2s))??
            warp_trl_mtx(itrl, 12) = warp_trl_mtx(itrl, 7)...
                + sum(time_series_event_spacing(4:5));  
            
            % original timing
            lo_time = trl_mtx_1_hold(itrl);
            hi_time = trl_mtx_1_hold(itrl) + trl_mtx(itrl, 12);
            
            % new timing
            warp_lo_time = trl_mtx_1_hold(itrl) + warp_trl_mtx(itrl, 1);
            warp_hi_time = trl_mtx_1_hold(itrl) + warp_trl_mtx(itrl, 12);

        else
            
            % load event times
            warp_trl_mtx(itrl, trl_mtx_col_seq(itsq)) = ...
                warp_trl_mtx(itrl, trl_mtx_col_seq(itsq-1)) + time_series_event_spacing(itsq-1);
            
            % original timing
            lo_time =  trl_mtx_1_hold(itrl) + trl_mtx(itrl, trl_mtx_col_seq(itsq-1));
            hi_time = trl_mtx_1_hold(itrl) + trl_mtx(itrl, trl_mtx_col_seq(itsq));
            
            % new timing
            %[trl_mtx_1_hold(itrl) trl_mtx_col_seq(itsq-1) warp_trl_mtx(itrl, trl_mtx_col_seq(itsq-1))]
            warp_lo_time =  trl_mtx_1_hold(itrl) + warp_trl_mtx(itrl, trl_mtx_col_seq(itsq-1));
            warp_hi_time = trl_mtx_1_hold(itrl) + warp_trl_mtx(itrl, trl_mtx_col_seq(itsq));

        end
        
        % original chunk frames
        chunkft_idx = trialft >= lo_time & trialft < hi_time;
        chunkft = trialft(chunkft_idx);
        %
        if max(trialft) < lo_time
            warning('frame times do not match behavior times')
            
            %[min(trialft) max(trialft)]
            %[lo_time hi_time]
            
            trial_duration  = (max(trialft) - (hi_time - lo_time));
            if trial_duration <=0
                trial_duration = 0.01;
            end
            %[max(trialft)-trial_duration max(trialft)]
            
            
            chunkft_idx = trialft >= max(trialft)-trial_duration & trialft < max(trialft);
            chunkft = trialft(chunkft_idx);
            
        elseif isempty(chunkft)

            chunkft = trialft(find(trialft >= lo_time, 1, 'first'));
            chunkft_idx = trialft==chunkft;
            if isempty(chunkft)
                chunkft = trialft(find(trialft <= hi_time, 1, 'last'));
                chunkft_idx = trialft==chunkft;
            end
        end
        %}
        
        % preprocess
        chunk_traces = trial_traces(:, chunkft_idx);
        warp_chunkft = linspace(warp_lo_time, warp_hi_time-0.01, nncf);
        
        % update frame times
        warp_frame_times_trl = [warp_frame_times_trl warp_chunkft];
        
  
        
        % iterate through each cell
        %chunkft
        %[min(chunkft), max(chunkft), length(warp_chunkft)]
        
        dummy_chunkft = linspace(min(chunkft), max(chunkft), length(warp_chunkft));
        warp_traces_chunk = nan(size(traces,1), nncf);
        for ic = 1:size(traces,1)
    
            % interp traces to new number of frames
            traces_ic_chunk = chunk_traces(ic, :);      
            if length(chunkft) == 1
                warp_traces_chunk(ic,:) = repmat(traces_ic_chunk, 1, length(dummy_chunkft));
            else
                warp_traces_chunk(ic,:) = interp1(chunkft, traces_ic_chunk, dummy_chunkft);
            end

        end
        
        % load chunk activity from all cells
        warp_traces_trl = [warp_traces_trl warp_traces_chunk];
        
        
        % POST trial times
        %{
        if trl_mtx_col_seq(itsq)==12 || (trl_mtx_col_seq(itsq)==11 && isnan(warp_trl_mtx(itrl, 11)))
            
            % find original frame times during event chunk
            lo_time = sum(trl_mtx(itrl, [1 12]));
            hi_time = max(trialft);
            postft_idx = trialft >= lo_time & trialft <= hi_time;
            
            % load post trial traces
            warp_traces_trl = [warp_traces_trl trial_traces(:, postft_idx)];
            
            % post trial frame times
            warp_frame_times_post = sum(warp_trl_mtx(itrl, [1 12])) + (trialft(postft_idx) - sum(trl_mtx(itrl, [1 12])));
            warp_frame_times_trl = [warp_frame_times_trl warp_frame_times_post];
            
        end
        %}
        
        
        % break if no reward
        if trl_mtx_col_seq(itsq)==11 && isnan(warp_trl_mtx(itrl, 11))
            break
        end
        
    end
    
    % update trial index
    warp_trl_idx_trl = repmat(itrl, length(warp_frame_times_trl), 1);
    
    
    % load output
    warp_traces = [warp_traces warp_traces_trl]; 
    warp_frame_times = [warp_frame_times warp_frame_times_trl]; 
    warp_trl_idx = [warp_trl_idx; warp_trl_idx_trl]; 
    
    % reset column 1
    warp_frame_times(warp_trl_idx==itrl) = warp_frame_times(warp_trl_idx==itrl) - warp_trl_mtx(itrl,1);
    warp_trl_mtx(itrl,[6 9 10 7 1 11 12]) = warp_trl_mtx(itrl,[6 9 10 7 1 11 12]) - warp_trl_mtx(itrl,1);
    warp_trl_mtx(itrl,1) = trl_mtx_1_hold(itrl) + warp_trl_mtx(itrl,1);

    
end







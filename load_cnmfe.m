function [traces_out, frame_times_out, trl_idx_out] = load_cnmfe(fp)
% load essential output from cnmfe

%E:\Projects\InProgress\GenProb\data\two_tone\train_lovar_imaging\651049m1\LED_gen14_probe\05

% traces
load([fp '\cnmfe_out.mat'], 'C_raw', 'rev_traces');
rev_traces = logical(rev_traces);
traces = C_raw(rev_traces,:);


% correct curve of traces
%
for itrace = 1:size(traces)
    traces(itrace,:) = traces_mean_correct(traces(itrace,:), 300);
end
%}

% frame times
load([fp '\pkl_frametimes.mat']);
frame_times = all_pkl_frames;

% frame times
load([fp '\pkl_trl_idx.mat'], 'trl_idx');

% catch glitch
if size(traces,2)<size(frame_times,1)
    fp
    frame_times = frame_times(1:size(traces,2));
    trl_idx = trl_idx(1:size(traces,2));
end


% soft preallocate
traces_out = [];
frame_times_out = [];
trl_idx_out = [];

% interp each trial to 10ms
for itrl = unique(trl_idx)'

    % observed and redrawn frames
    frame_times_obs = frame_times(trl_idx==itrl);
    frame_times_new = floor(frame_times_obs(1)*100)/100 : 0.01 : ceil(frame_times_obs(end)*100)/100;
    traces_obs = traces(:,trl_idx==itrl);
    traces_new = nan(size(traces,1), length(frame_times_new));

    % frame times out
    frame_times_out = [frame_times_out frame_times_new];
    
    % interp
    for itrace = 1:size(traces_obs,1)
        traces_new(itrace,:) = interp1(frame_times_obs, traces_obs(itrace,:), frame_times_new);
    end

    % traces out
    traces_out = [traces_out traces_new];
    
    % trial index out
    trl_idx_out = [trl_idx_out; repmat(itrl, length(frame_times_new), 1)];
    
end

% remove nans
nnan_idx = ~isnan(traces_out(1,:));
frame_times_out = frame_times_out(nnan_idx);
traces_out = traces_out(:,nnan_idx);
trl_idx_out = trl_idx_out(nnan_idx);










function ALL_image_mean_activity_timewarp(file_keywords, first_last)
% combine cells from a single stage from all animals

% fixed trial timeline: nose_poke onset, nose_poke offset, head_entry, tone
% on, reward, trial end
time_series_event_spacing = [0.5 1.1 1.5 2.0 3.5 2.0];

% get all files
fp = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_imaging\';
session_files = get_file_paths_targeted(fp, file_keywords);

% find unique subjects
subjs = [];
for isf = 1:size(session_files,1)
    subjs = [subjs; {session_files{isf}(strfind(session_files{isf}, '\\')+2 : strfind(session_files{isf}, '\\')+9)}];
end
subjs = unique(subjs);

% find one file per subject
session_files_hold = [];
for isubj = 1:length(subjs)
    subj_sf = session_files(contains(session_files, subjs{isubj}));
    if strcmp(first_last, 'first')
        session_files_hold = [session_files_hold; subj_sf(1)];
    elseif strcmp(first_last, 'last')
        session_files_hold = [session_files_hold; subj_sf(end)];
    elseif strcmp(first_last, 'all')
        session_files_hold = [session_files_hold; subj_sf(:)];
    else
        error('do better with first_last input')
    end
end
session_files = session_files_hold


% iterate through session files
all_mtx = [];
all_sort = [];
for isf = 1:size(session_files,1)

    % load session file
    session_files{isf}
    load(session_files{isf}, 'frame_times', 'traces', 'trl_idx', 'trl_mtx')
    
    if ~exist('trl_idx', 'var') || ~exist('traces', 'var')
        continue
    end
    
    % trial type
    tt = 'rwd'
    if strcmp(tt, 'probe')
        trial_nums = intersect(unique(trl_idx), find(trl_mtx(:,3)==0));
    elseif strcmp(tt, 'no rwd')
        trial_nums = intersect(unique(trl_idx), find(isnan(trl_mtx(:,11))));
    elseif strcmp(tt, 'rwd')
        trial_nums = intersect(unique(trl_idx), find(~isnan(trl_mtx(:,11))));
    end
    
    % compute matrices and order    
    [sf_mtx,sf_srt_idx] = image_mean_activity_timewarp...
        (trl_mtx, trl_idx, frame_times, traces, 1:size(traces,1),...
        trial_nums);

    % close figures
    close; close
    
    % load
    all_mtx = [all_mtx; sf_mtx];
    all_sort = [all_sort; sf_srt_idx];
    
    % clear lingering variables
    clearvars('frame_times', 'traces', 'trl_idx', 'trl_mtx')
end

% sort rows by max pixel
max_pix = nan(size(all_mtx,1),1); 
for irow = 1:size(all_mtx,1)
    max_pix(irow) = find(all_mtx(irow,:)==max(all_mtx(irow,:))); 
end
[~, sort_idx] = sort(max_pix);
all_mtx = all_mtx(sort_idx,:);

% plot
figure; hold on
imagesc(all_mtx)

% red line
event_frame = cumsum([0 time_series_event_spacing]).*100;

xticks_hold = [];
for i = event_frame 
   plot(i.*[1 1], ylim, 'r-', 'linewidth', 1.5) 
   xticks_hold = [xticks_hold i];
end

% aesthetics
set(gca, 'YDir','reverse')
ylim([0+0.5 size(all_mtx,1)+0.5])
xlim([0.5 size(all_mtx,2)+0.5])
xticks(xticks_hold)
xticklabels_universe = {'NP1', 'NP2', 'HE1', 'TO', 'RD', 'RWD', 'TE'};
xticklabels(xticklabels_universe)

% lineplot
figure; hold on 
plot(smooth(mean(all_mtx), 10), 'k-', 'linewidth', 2)
ylim([0 1])
set(gca,'TickLength',[0, 0]); box off;
%histogram(max_pix, 0 : 25 : size(all_mtx,2)+1, 'Normalization', 'probability', 'DisplayStyle', 'stairs')

% red line
event_frame = cumsum([0 time_series_event_spacing]).*100;
for i = event_frame 
   plot(i.*[1 1], ylim, 'r-', 'linewidth', 1.5) 
end
xlim([0.5 size(all_mtx,2)+0.5])




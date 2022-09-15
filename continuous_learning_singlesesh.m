function [wait_prefs, rich_waits, poor_waits] = continuous_learning_singlesesh(session_path, wdw_size)
% plot accuracy over a sliding window that runs seamlessly across training
% sessions

% extract subject id
slash_idx = strfind(session_path, '\');
subject_id = session_path(slash_idx(end-1)+1 : slash_idx(end)-1);

% training group string
subj_pos = strfind(session_path, subject_id);
slash_pos = strfind(session_path, '\');
slash_pos = slash_pos(slash_pos<subj_pos);
slash_pos = slash_pos(end-1: end);
train_group = session_path(slash_pos(1)+1 : slash_pos(2)-1);
train_group(strfind(train_group, '_')) = ' ';

%% combine trl_mtxs
[all_trl_mtx, session_number, problem_number] = ALL_trl_mtx(session_path);

% only probe trials
session_number_hold = session_number;
problem_number_hold = problem_number;
session_number = session_number(all_trl_mtx(:,3)==0);
problem_number = problem_number(all_trl_mtx(:,3)==0);
all_trl_mtx_hold = all_trl_mtx;
all_trl_mtx = all_trl_mtx(all_trl_mtx(:,3)==0,:);

% all rich tones
load('unqfrq41.mat', 'unqfrq41')
if contains(session_path, 'mevar')
    rich_tones = unqfrq41(13:29);
elseif contains(session_path, 'hivar')
    rich_tones = unqfrq41([5 13 21 29 37]);
end
%rich_tones = unqfrq41([5 13 16 18 19 20 21 22 23 24 26 29 37]);
rich_tone_idx = ismember(floor(all_trl_mtx_hold(:,2)), rich_tones);
probe_idx = all_trl_mtx_hold(:,3)==0;


%% sliding window
rwt_last = all_trl_mtx_hold(rich_tone_idx & all_trl_mtx_hold(:,3)==0, 12);
    rwt_last = rwt_last(1);
pwt_last = all_trl_mtx_hold(~rich_tone_idx & all_trl_mtx_hold(:,3)==0, 12); 
    pwt_last = pwt_last(1);
wait_prefs = nan(size(all_trl_mtx_hold,1),1);
rich_waits = nan(size(all_trl_mtx_hold,1),1);
poor_waits = nan(size(all_trl_mtx_hold,1),1);

for islide = 1:size(all_trl_mtx_hold,1)
    
    win_lo = islide-floor(wdw_size/2); if win_lo<1; win_lo=1; end
    win_hi = islide+floor(wdw_size/2); if win_hi>size(all_trl_mtx_hold,1); win_hi=size(all_trl_mtx_hold,1); end
    
    islide_idx = ismember(1:size(all_trl_mtx_hold,1), win_lo:win_hi);

    if length(unique(session_number_hold(islide_idx)))>1
        wait_prefs(islide+round(wdw_size/2)) = nan;
        continue
    end
    
% rich wait times
rwt = mean(all_trl_mtx_hold(islide_idx' & rich_tone_idx & probe_idx, 12));
if isnan(rwt)
    rwt = rwt_last;
end

% poor wait times
pwt = mean(all_trl_mtx_hold(islide_idx' & ~rich_tone_idx & probe_idx, 12));
if isnan(pwt)
    pwt = pwt_last;
end

% LOAD proportion (rich-poor / rich+poor)
wait_prefs(islide) = (rwt-pwt) / (rwt+pwt);
rich_waits(islide) = rwt;
poor_waits(islide) = pwt;

% reset lasts
rwt_last = rwt;
pwt_last = pwt;

end

end
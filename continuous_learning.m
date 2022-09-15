function continuous_learning(wdw_size, subject_id)
% plot accuracy over a sliding window that runs seamlessly across training
% sessions

%get_file_paths_targeted('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_consistent', 'train', '658103m3')

% grey and red colors
colors_gr = [.6 .6 .6; [141 2 31]./255];
colors_diff = diff(colors_gr);

%% all paths
session_paths = [get_file_paths_targeted(...
    'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data', 'novar0', subject_id);...
    get_file_paths_targeted('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data', 'lovar0', subject_id);...
    get_file_paths_targeted('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data', 'mevar0', subject_id);...
    get_file_paths_targeted('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data', 'hivar0', subject_id);...
    get_file_paths_targeted('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data', 'ctl0', subject_id);...
    get_file_paths_targeted('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data', 'exvar0', subject_id)];
%session_paths = session_paths(end)
%session_paths = session_paths(contains(session_paths,'meta'));
session_paths = session_paths(~contains(session_paths, 'covid_hold'));

% training group string
fp = session_paths{1};
subj_pos = strfind(fp, subject_id);
slash_pos = strfind(fp, '\');
slash_pos = slash_pos(slash_pos<subj_pos);
slash_pos = slash_pos(end-1: end);
train_group = fp(slash_pos(1)+1 : slash_pos(2)-1);
session_paths = session_paths(contains(session_paths, train_group));
train_group(strfind(train_group, '_')) = ' ';

%% combine trl_mtxs
[all_trl_mtx, session_number, problem_number] = ALL_trl_mtx(session_paths);

%{
delete_prob_nums = [5 max(session_number(problem_number==5));...
    6 max(session_number(problem_number==6))];
for idpn = 1:size(delete_prob_nums,1)
    delete_idx = problem_number==delete_prob_nums(idpn,1) ...
        & session_number==delete_prob_nums(idpn,2);
    
    session_number(delete_idx,:) = [];
    problem_number(delete_idx,:) = [];
    all_trl_mtx(delete_idx,:) = [];

end
%}


% only probe trials
session_number_hold = session_number;
problem_number_hold = problem_number;
session_number = session_number(all_trl_mtx(:,3)==0);
problem_number = problem_number(all_trl_mtx(:,3)==0);
all_trl_mtx_hold = all_trl_mtx;
all_trl_mtx = all_trl_mtx(all_trl_mtx(:,3)==0,:);
all_trl_mtx_hold_row= 1:size(all_trl_mtx_hold,1);



% add empty sessions between each problem??


% problem number colors
%pn_colors = summer(6);
pn_colors = bone(10); pn_colors = pn_colors(2:7,:);

% all rich tones
load('unqfrq41.mat', 'unqfrq41')
%rich_tones = unqfrq41(13:29);
rich_tones = unqfrq41([5 13 16 18 19 20 21 22 23 24 26 29 37]);
rich_tone_idx = ismember(floor(all_trl_mtx_hold(:,2)), rich_tones);
probe_idx = all_trl_mtx_hold(:,3)==0;


%% sliding window
rwt_last = all_trl_mtx_hold(rich_tone_idx & all_trl_mtx_hold(:,3)==0, 12); rwt_last = rwt_last(~isnan(rwt_last)); rwt_last = rwt_last(1);
pwt_last = all_trl_mtx_hold(~rich_tone_idx & all_trl_mtx_hold(:,3)==0, 12); pwt_last = pwt_last(~isnan(pwt_last)); pwt_last = pwt_last(1);
wait_prefs = nan(size(all_trl_mtx_hold,1),1);


    %figure; plot(all_trl_mtx_hold_row(rich_tone_idx & all_trl_mtx_hold(:,3)==0), all_trl_mtx_hold(rich_tone_idx & all_trl_mtx_hold(:,3)==0, 12), 'o')
    %hold on; plot(all_trl_mtx_hold_row(~rich_tone_idx & all_trl_mtx_hold(:,3)==0), all_trl_mtx_hold(~rich_tone_idx & all_trl_mtx_hold(:,3)==0, 12), 'o')

    
for islide = 1:size(all_trl_mtx_hold,1)
    
    win_lo = islide-floor(wdw_size/2); if win_lo<1; win_lo=1; end
    win_hi = islide+floor(wdw_size/2); if win_hi>size(all_trl_mtx_hold,1); win_hi=size(all_trl_mtx_hold,1); end
    
    islide_idx = ismember(1:size(all_trl_mtx_hold,1), win_lo:win_hi);

    if length(unique(session_number_hold(islide_idx)))>1
        wait_prefs(islide+round(wdw_size/2)) = nan;
        continue
    end
    
    % rich wait times
    rwt = nanmean(all_trl_mtx_hold(islide_idx' & rich_tone_idx & probe_idx, 12));
    if isnan(rwt)
        rwt = rwt_last;
    end

    % poor wait times
    pwt = nanmean(all_trl_mtx_hold(islide_idx' & ~rich_tone_idx & probe_idx, 12));
    if isnan(pwt)
        pwt = pwt_last;
    end


    % proportion (rich-poor / rich+poor)
    wait_prefs(islide) = (rwt-pwt) / (rwt+pwt);

    % reset lasts
    rwt_last = rwt;
    pwt_last = pwt;

end

%% plot

% smooth and interp wait prefs
wait_prefs_int = nanfastsmooth(wait_prefs, 20, 1, 1);

% figure
hold on

% plot problem background colors
for ipn = unique(problem_number)'
   
    pn_start = find(problem_number_hold==ipn, 1, 'first');
    pn_end = find(problem_number_hold==ipn, 1, 'last');
    
    %pn_color_local = pn_colors(ipn,:);
    pn_color_local = [1 1 1];
    rectangle('Position', [pn_start -1 pn_end-pn_start 2], 'FaceColor', pn_color_local, 'EdgeColor', [0 0 0], 'LineWidth', 2)
end

% plot curve and dots
%plot(1:length(wait_prefs), wait_prefs, 'o', 'color', .7.*[1 1 1])
cutoff = -0.1;
for idot = 1:length(wait_prefs)
    if isnan(wait_prefs(idot)); continue; end
    %
    if wait_prefs(idot)<cutoff 
        current_color = colors_gr(1,:);
    else
        current_color = colors_gr(1,:)+colors_diff.*((wait_prefs(idot)-cutoff)/(1-cutoff));
    end
    %}
    %current_color = pn_colors(problem_number_hold(idot),:);
    plot(idot, wait_prefs(idot), 'o', 'color', current_color)
end
plot(1:length(wait_prefs_int), wait_prefs_int, '-', 'color', 0.25.*[1 1 1], 'linewidth', 3)


% session bounds
ylim([-1 1])
[sesh_bnds, ~] = hist(session_number_hold, unique(session_number_hold));
sesh_bnds = cumsum(sesh_bnds);
plot([sesh_bnds;sesh_bnds]+0.5, ylim, 'k-')

% chance
plot(xlim, [1 1].*0, 'k--')

% aesthetics
xlim([1 length(wait_prefs)])
title([subject_id '; ' train_group])
ylabel('Preference for rich tone')
xticks([])
xlabel('Problem, Session')
set(gca,'TickLength',[0, 0])

end
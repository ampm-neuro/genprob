function [within_sesh_zchange, between_sesh_zchange, all_within_diffs, all_between_diffs] = wthn_btwn_learn(varargin)
% compares learning between session halves (within-session learning) with 
% learning between sessions (between-session learning)
%
% input cell of string(s) that can be used to distinguish between file 
% names

% umbrella path
folderpath = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\Gen_richards';

% constrain via string inputs
fpaths = path_comp(folderpath, varargin);

%unique subjects-stage ids
folderpath_slash = strfind(folderpath, '\');
fpath_slash = strfind(fpaths, '\');
unq_subjs = cell(length(fpaths),1);
unq_stage = cell(length(fpaths),1);
comb_subj_stage = cell(length(fpaths),1);
for ifpath = 1:length(fpath_slash)
    bound_slash = find(fpath_slash{ifpath}>max(folderpath_slash),2, 'first');
    fpath_subj_bounds = fpath_slash{ifpath}(bound_slash) + [1 -1];
    unq_subjs{ifpath} = fpaths{ifpath}(fpath_subj_bounds(1):fpath_subj_bounds(2));
    unq_stage{ifpath} = fpaths{ifpath}(fpath_subj_bounds(2)+2 : fpath_subj_bounds(2)+6);
    comb_subj_stage{ifpath} = [unq_subjs{ifpath} unq_stage{ifpath}];
end
[~,~,subj_stage_id] = unique(comb_subj_stage);


% divide sessions into halves
session_halves = cell(length(fpaths),2);
for isesh = 1:length(fpaths)  
    
    %load trial matrix
    load(fpaths{isesh}, 'trl_mtx', 'medass_cell');
    
    %window size
    num_trls = size(trl_mtx,1);
    winsz = ceil(num_trls/2);
    
    %preallocate
    rich_waits = cell(1);
    poor_waits = cell(1);

    %for each window
    window_count = 0;
    for iw = 1:winsz:size(trl_mtx,1)
        
        window_count = window_count+1;
        
        %trl mtx for this sesh and window
        if num_trls>=iw+winsz
            first_bin = iw;
            last_bin = iw+winsz-1;
        else
            first_bin = iw;
            last_bin = num_trls;
        end
       
        local_trl_mtx = trl_mtx(first_bin:last_bin,:);
        
        %compute wait durations
        [~, wait_durations, ~, ~, p_dist] = wait_times(local_trl_mtx, medass_cell, 0);
        
        min_samps = 1;
        rich_waits = cell2mat(wait_durations(p_dist>0.50));       
            if length(rich_waits) < min_samps
                rich_waits = nan; 
            end
            
        poor_waits = cell2mat(wait_durations(p_dist<0.50));
            if length(poor_waits) < min_samps
                poor_waits = nan; 
            end
               
        %load rich and poor for each half
        session_halves{isesh, window_count} = [{rich_waits} {poor_waits}];
            
        %break    
        if last_bin == size(trl_mtx,1)
            break
        end
        
    end
end

% session half zdiffs
session_half_zdiffs = nan(size(session_halves,1),2);
for isesh = 1:size(session_halves,1)
    %first half
    if length(session_halves{isesh,1}{1})<min_samps
        session_half_zdiffs(isesh,1) = nan;
    else
        session_half_zdiffs(isesh,1) = zdiff(session_halves{isesh,1}{1}, session_halves{isesh,1}{2});
    end

    %second half
    if length(session_halves{isesh,1}{2})<min_samps
        session_half_zdiffs(isesh,2) = nan;
    else
        session_half_zdiffs(isesh,2) = zdiff(session_halves{isesh,2}{1}, session_halves{isesh,2}{2});
    end
end

% zdiffs
within_sesh_zchange = [];
    all_within_diffs = [];
between_sesh_zchange = [];
    all_between_diffs = [];
for isubstg = unique(subj_stage_id)'
    local_zdiffs = session_half_zdiffs(subj_stage_id==isubstg,:);

    % within-session zdiff change
    within_diff = local_zdiffs(:,2) - local_zdiffs(:,1);
    within_diff = within_diff';
    [within_diff, all_within_diffs] = match_mtx_width(within_diff, all_within_diffs);
    all_within_diffs = [all_within_diffs; within_diff];
    sum_wthn_sesh = nansum(within_diff);
    within_sesh_zchange = [within_sesh_zchange; sum_wthn_sesh];

    %between-session zdiff change
    third_half = local_zdiffs(2:end,1);
    between_diff = third_half - local_zdiffs(1:end-1,2);
    between_diff = between_diff';
    [between_diff, all_between_diffs] = match_mtx_width(between_diff, all_between_diffs);
    all_between_diffs = [all_between_diffs; between_diff];

    if sum(~isnan(between_diff))>0
        sum_btwn_sesh = nansum(between_diff);
    elseif size(local_zdiffs,1)==1            
        sum_btwn_sesh = 0;
    else
        sum_btwn_sesh = nan;
    end
    sum_btwn_sesh = sum_btwn_sesh;
    between_sesh_zchange = [between_sesh_zchange; sum_btwn_sesh];
    
    
    %sum_total_learn = nansum([sum_wthn_sesh sum_btwn_sesh])
    
end


%CUMULATIVE
%{
all_within_diffs_cum = nancumsum(all_within_diffs,2);
    all_within_diffs_cum(isnan(all_within_diffs)) = nan;
    all_within_diffs = all_within_diffs_cum;
all_between_diffs_cum = nancumsum(all_between_diffs,2);
    all_between_diffs_cum(isnan(all_between_diffs)) = nan;
    all_between_diffs = all_between_diffs_cum;
%}

%plot learning curve across single stage
%{
colors = distinguishable_colors(2);
figure; hold on
errorbar(1:size(all_within_diffs,2), nanmean(all_within_diffs,1), nanstd(all_within_diffs,[],1)./sqrt(sum(~isnan(all_within_diffs),1)), '-o', 'color', colors(1,:));
errorbar(1.5:1:size(all_within_diffs,2)-0.5, nanmean(all_between_diffs,1), nanstd(all_between_diffs,[],1)./sqrt(sum(~isnan(all_between_diffs),1)), '-o', 'color', colors(2,:));
xlim([.5 14.5])
set(gca,'TickLength',[0, 0]); box off;
hold on; plot(xlim, [1 1].*0, 'k--')
legend({'Within-sesssion learning', 'Between-session learning'}, 'location', 'northeastoutside')
ylim([-6.5 6.55])
xticks(1:14)
%}

% plot errorbar
%{
figure; hold on
bar([nanmean(within_sesh_zchange) nanmean(between_sesh_zchange)])
errorbar_scatter([{within_sesh_zchange} {between_sesh_zchange}]);
errorbar([nanmean(within_sesh_zchange) nanmean(between_sesh_zchange)], ...
    [nanstd(within_sesh_zchange) nanstd(between_sesh_zchange)]...
    ./sqrt([sum(~isnan(within_sesh_zchange)) sum(~isnan(between_sesh_zchange))]), 'k.')
set(gca,'TickLength',[0, 0]); box off;

[~,p_within_sesh,~,stats_within_sesh] = ttest(within_sesh_zchange);
sig_asterisks(p_within_sesh, 1, max([nanmean(within_sesh_zchange) nanmean(between_sesh_zchange)])+0.5)
[~,p_between_sesh,~,stats_between_sesh] = ttest(between_sesh_zchange);
sig_asterisks(p_between_sesh, 2, max([nanmean(within_sesh_zchange) nanmean(between_sesh_zchange)])+0.5)
[~,p_ind,~,stats_ind] = ttest2(within_sesh_zchange, between_sesh_zchange);
sig_asterisks(p_ind, 1.5, max([nanmean(within_sesh_zchange) nanmean(between_sesh_zchange)])+0.75)

%ylim_hold = ylim;
%ylim([-1 ylim_hold(2)+1])
xticks([1 2])
xticklabels({'WithinSesh', 'BetweenSesh'})
ylabel('Preference for rich tone (zdiff)')
%}

end

function [mtx1, mtx2] = match_mtx_width(mtx1, mtx2)
    
    if isempty(mtx1) || isempty(mtx2)
        return
    end
    
    width1 = size(mtx1,2);
    width2 = size(mtx2,2);
    
    if width1 < width2
        mtx1 = [mtx1 nan(size(mtx1,1), width2-width1)];
    elseif size(mtx1,2) > size(mtx2,2)
        mtx2 = [mtx2 nan(size(mtx2,1), width1-width2)];
    end
end

function B=nancumsum(A,dim,nmode)
% NANCUMSUM: Cumulative sum of a matrix, with user-specified treatment of NaNs.
%     Computes the cumulative sum of matrix A along dimension DIM, allowing
%     the user to replace NaNs with zeros, to skip over them, or to reset
%     on NaNs, maintaining NaNs as placeholders. 
% 
% USAGE: B = nancumsum(A, DIM, NMODE)
%
% ARGUMENTS:
%
% A:    Input matrix.
%
% B:    Output cumulative sum matrix, treating NaNs as determined by nmode.
%
% DIM:  B = nancumsum(A, DIM) returns the nan-cumulative sum of the elements
%       along the dimension of A specified by scalar DIM. For example,
%       nancumsum(A,1) works down the columns, nancumsum(A,2) works
%       across the rows. If DIM is not specified, it defaults to the first
%       non-singleton dimension of A. 
%
% NMODE: specifies how NaNs should be treated. Acceptable values are:
%       1: REPLACE NaNs with zeros (default).
%       2: MAINTAIN NaNs as position holders in B. (Skip NaNs without reset.)
%       3: RESET sum on NaNs, replacing NaNs with zeros.
%       4: RESET sum on NaNs, maintaining NaNs as position holders.
%
% EXAMPLES:
%
% 1) a = [NaN,2:5];
%
% nancumsum(a)
% ans =
%     0     2     5     9    14
%
% nancumsum(a,[],2)
% ans =
%   NaN     2     5     9    14
%
% nancumsum(a,[],3)
% ans =
%     2     5     9    14
%
% 2) a = magic(3); a(5)=NaN;
%
% b = nancumsum(a,2) % (Default NMode = 1)
% b =
%     8     9    15
%     3     3    10
%     4    13    15
%
% b = nancumsum(a,2,2)
% b =
%     8     9    15
%     3   NaN    10
%     4    13    15
%
% b = nancumsum(a,2,3)
% b =
%     8     9    15
%     3     0     7
%     4    13    15
%
% b = nancumsum(a,2,4)
% b =
%     8     9    15
%     3   NaN     7
%     4    13    15
% See also: cumsum, nansum, nancumprod, nanmean, nanmedian, ...
% (nancumprod is available from the FEX. Other nan* may require Toolboxes)
% Brett Shoelson
% brett.shoelson@mathworks.com
% 05/04/07
%
% Revision: 08/28/11
% Fixed bug in option 2 (faulty reset). Thanks to Andrew Stevens and Rick
% Patterson for reporting it. Also, eliminated old option 3 (deleting NaNs
% in a vector) as a trivial case and added two new options. 
% Revision: 09/15/11
% Fixed a bug with multiple NaNs. Thanks to Tim Yates.
%
% Copyright The MathWorks, Inc. 2011
% Set defaults, check and validate inputs
if nargin < 3
    nmode = 1;
end
if ~ismember(nmode,1:4)
    error('NANCUMSUM: unacceptable value for nmode parameter.');
end
if nargin < 2 || isempty(dim)
    if ~isscalar(A)
        dim = find(size(A)>1);
        dim = dim(1);
    else
        % For scalar inputs (no nonsingleton dimension)
        dim = 1;
    end
end
% Calculate cumulative sum, depending on selection of nmode
switch nmode
    case 1
        % TREAT NaNs as 0's
        B = A;
        B(B~=B) = 0;
        B = cumsum(B, dim);
    case 2
        % DO NOT INCREMENT, BUT USE NaNs AS PLACEHOLDERS.
        B = nancumsum(A, dim, 1);
        B(A~=A) = NaN;
     case 3
        % RESET sum on NaNs, replacing NaNs with zeros.
        naninds = find(A~=A);
        for ii = 1:numel(naninds)
            B = nancumsum(A, dim, 1);
            A(naninds(ii)) = -B(naninds(ii));
        end
        B = cumsum(A,dim);
    otherwise %case 4
        % RESET sum on NaNs, maintaining NaNs as position holders.
        naninds = find(A~=A);
        B = nancumsum(A,dim,3);
        B(naninds)= NaN;
end

end

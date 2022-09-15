function [discrim_optoOFF, discrim_optoON, session_file] = opto_discrim_singlesubj(subject, problem, first_last, varargin)
% find discrimination index (cohen D) for rich an poor tones in a single
% subject under optoON and optoOFF conditions


%% input file folder
if nargin == 4
    fp = varargin{1};
else
    fp = ['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\'];
end



%% get session

% opto session files from this subject & problem
session_files = get_file_paths_targeted(fp, {subject, ['var0' num2str(problem) '_'], '_opto-'})

% if empty, return
if isempty(session_files)
    discrim_optoOFF = nan;
    discrim_optoON = nan;
    session_file = 'no file';
    return
end

% first last
if first_last == 1
    session_file = session_files{1};
elseif first_last == 2 && size(session_files,1)>=2
    session_file = session_files{end};
else
    warning('no valid file')
end



%% compute discrimination

% load session file
load(session_file, 'trl_mtx')
            
% idx
norwd_idx = trl_mtx(:,3)==0;
optoON_idx = trl_mtx(:,13)==1;
rich_idx = rich_trl_idx(trl_mtx); 

% minimum rich data points
min_rich_norwd_trials = 1;
if length(trl_mtx(norwd_idx & optoON_idx & rich_idx, 12)) < min_rich_norwd_trials
    discrim_optoON = nan;
end

% tones
poor_optoOFF = trl_mtx(norwd_idx & ~optoON_idx & ~rich_idx, 12);
rich_optoOFF = trl_mtx(norwd_idx & ~optoON_idx & rich_idx, 12);
poor_optoON = trl_mtx(norwd_idx & optoON_idx & ~rich_idx, 12);
rich_optoON = trl_mtx(norwd_idx & optoON_idx & rich_idx, 12);

% cohen d (optoOFF, optoON)
discrim_optoOFF = computeCohen_d(rich_optoOFF, poor_optoOFF);
discrim_optoON = computeCohen_d(rich_optoON, poor_optoON);




function [all_trl_mtx, session_number, problem_number] = ALL_trl_mtx(session_paths)
% combined trl_mtx s from multiple sessions
% session_number provides a session id for each row

% preallcoate
all_trl_mtx = [];
session_number = [];
problem_number = [];

%iterate through sessions
for isesh = 1:size(session_paths,1)
    
    % path string
    if iscell(session_paths)
        path_string = session_paths{isesh};
    else
        path_string = session_paths;
    end
    %load data file
    load(path_string, 'trl_mtx')

    % zscore probe wait times
    %trl_mtx(trl_mtx(:,3)==0,12) = zscore_mtx(trl_mtx(trl_mtx(:,3)==0,12))-2;
    
    %load output
    %size(all_trl_mtx)
    %size(trl_mtx)
    
    if size(trl_mtx,2) < size(all_trl_mtx,2)
        trl_mtx = [trl_mtx nan(size(trl_mtx,1), size(all_trl_mtx,2)-size(trl_mtx,2))];
    elseif size(trl_mtx,2) > size(all_trl_mtx,2)
        all_trl_mtx = [all_trl_mtx nan(size(all_trl_mtx,1), size(trl_mtx,2)-size(all_trl_mtx,2))];
    end
    
    all_trl_mtx = [all_trl_mtx; trl_mtx];
    session_number = [session_number; repmat(isesh, size(trl_mtx,1), 1)];
    
    if contains(path_string, 'novar0')
        pn = str2num(path_string(strfind(path_string, 'novar0')+length('novar0')));
    elseif contains(path_string, 'lovar0')
        pn = str2num(path_string(strfind(path_string, 'lovar0')+length('lovar0')));
    elseif contains(path_string, 'mevar0')
        pn = str2num(path_string(strfind(path_string, 'mevar0')+length('mevar0')));
    elseif contains(path_string, 'hivar0')
        pn = str2num(path_string(strfind(path_string, 'hivar0')+length('hivar0')));
    elseif contains(path_string, 'ctl0')
        pn = str2num(path_string(strfind(path_string, 'ctl0')+length('ctl0')));
    elseif contains(path_string, 'exvar0')
        pn = str2num(path_string(strfind(path_string, 'exvar0')+length('exvar0')));
    elseif contains(path_string, 'probe_0')
        pn = str2num(path_string(strfind(path_string, 'probe_0')+length('probe_0')));
    elseif contains(path_string, 'probe_opto_0')
        pn = str2num(path_string(strfind(path_string, 'probe_opto_0')+length('probe_opto_0')));
    else
        error('What learning condition is this?')
    end
    problem_number = [problem_number; repmat(pn, size(trl_mtx,1), 1)];

end


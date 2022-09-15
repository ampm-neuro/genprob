function [topDown_idx, bottomUp_idx] = cell_type_idx(cell_regist_mtx)
% outputs logical matrices with dimensions of cell_regist_mtx that indicate
% whether that cell is a topDown cell or a bottomUp cell

% session types
probe_sessions = [1:3:19 20];
problem_sessions = [2:3:17; 3:3:18];

%% top down index
tdi = zeros(size(cell_regist_mtx,1), 6);

% iterate through cells
for icell = 1:size(cell_regist_mtx,1)
    % iterate through probes
    for iprobe = 1:length(probe_sessions)-2
        
        % corresponding problem sessions
        corresp_problem_sessions = probe_sessions(iprobe)+1:probe_sessions(iprobe)+2;
        
        % if active on this probe, and subsequent problem, and final probe8
        if cell_regist_mtx(icell, probe_sessions(iprobe))>0 && sum(cell_regist_mtx(icell,corresp_problem_sessions),2)>0 && cell_regist_mtx(icell, probe_sessions(end))>0
            tdi(icell, iprobe) = 2;
        elseif cell_regist_mtx(icell, probe_sessions(iprobe))>0 && sum(cell_regist_mtx(icell,corresp_problem_sessions),2)>0
            tdi(icell, iprobe) = 1;
        end
        
    end
end

% load
topDown_idx = zeros(size(cell_regist_mtx));
topDown_idx(:,1:3:16) = tdi;



%% bottom up index
bui = zeros(size(cell_regist_mtx));

% iterate through cells
for icell = 1:size(cell_regist_mtx,1)
    
    origin_sesh = find(cell_regist_mtx(icell,:)>0,1,'first');
    
    % iterate through problem numbers
    for iproblem = 1:size(problem_sessions,2)

        subs_probe = min(probe_sessions(probe_sessions>max(problem_sessions(:,iproblem))));
        
        % iterate through each session number in this problem
        for session_number = problem_sessions(:,iproblem)'
            
            % if originated on this problem, active on subsequent probe, and final probe8
            if origin_sesh==session_number && cell_regist_mtx(icell,subs_probe)>0 && cell_regist_mtx(icell, probe_sessions(end))>0
                bui(icell, session_number) = 4;
                
            % if originated on this problem and active on subsequent probe
            elseif origin_sesh==session_number && cell_regist_mtx(icell,subs_probe)>0
                bui(icell, session_number) = 3;
                
            % if originated on subsequent probe, and active in final probe8
            elseif cell_regist_mtx(icell,subs_probe)>0 && cell_regist_mtx(icell, probe_sessions(end))>0
                %bui(icell, session_number) = 2;
                
            % if originated on subsequent probe
            elseif cell_regist_mtx(icell,subs_probe)>0
                %bui(icell, session_number) = 1;
                
            end
        end
    end
end

% load
bottomUp_idx = bui;




function all_footprints_bottomUP_colored(footprint_mtx, cell_regist_mtx)
% plots footprints for all probes colored according to which probe
%
% try aligned_data_struct.spatial_footprints_corrected for spatial
% footprints
%

figure; hold on

% colors
%{
all_bcolors = [.3.*[255 255 255];...
                186 41 41; ...
                232 116 0; ...
                232 204 37; ...
                28 142 22; ...
                22 178 247; ...
                48 83 193]./255;
%}
            
all_bcolors = hsv(20);

% session order
session_order = [1 9 15 2 10 16 3 11 17 4 12 18 5 13 19 6 14 20 7 8];
cell_regist_mtx = cell_regist_mtx(:,session_order);
footprint_mtx = footprint_mtx(session_order);
probe_sessions = [1:3:19 20];
problem_sessions = setdiff(1:20, probe_sessions);

% session plotting coordinates
yspace = 275;
spc = ...
[0 0; 200 -1*yspace; 200 -2.2*yspace; ...
400 0; 600 -1*yspace; 600 -2.2*yspace; ...
800 0; 1000 -1*yspace; 1000 -2.2*yspace; ...
1200 0; 1400 -1*yspace; 1400 -2.2*yspace; ...
1600 0; 1800 -1*yspace; 1800 -2.2*yspace; ...
2000 0; 2200 -1*yspace; 2200 -2.2*yspace; ...
2400 0; 2800 0];




%% plot
%
for icell = 1:size(cell_regist_mtx,1)
    
    % first session for this cell (term sessions will only plot if cell
    % active)
    origin_sesh = find(cell_regist_mtx(icell,:)>0, 1, 'first');
    terminal_probe_sesh = min(probe_sessions(probe_sessions>=origin_sesh));
    terminal_probe_num = find(probe_sessions==min(probe_sessions(probe_sessions>=origin_sesh)),1);
    active_during_terminal = cell_regist_mtx(icell, terminal_probe_sesh)>0;
    
    % if origin session is a problem session, color from origin to term
    % probe (if active)
    %if ismember(origin_sesh, problem_sessions)
    if origin_sesh>1
        color_rng = origin_sesh:terminal_probe_sesh;
    else
        color_rng = [];
    end
    
    

    for isesh = 1:size(cell_regist_mtx,2)-1

        current_cell = cell_regist_mtx(icell,isesh);

        if current_cell==0
            continue
        end

        if ismember(isesh, color_rng) && origin_sesh~=terminal_probe_sesh && active_during_terminal == 1
            %plot_color = all_bcolors(terminal_probe_num, :);
            plot_color = all_bcolors(origin_sesh, :);
        elseif ismember(isesh, color_rng) && origin_sesh==terminal_probe_sesh && origin_sesh>1
            plot_color = all_bcolors(origin_sesh, :); %probe originating
        else
            plot_color = .8.*[1 1 1];
        end
        
        linewidth = 1;
        %{
        if ismember(isesh, color_rng) && origin_sesh~=terminal_probe_sesh && active_during_terminal == 1
            linewidth = 2;
        else
            linewidth = 1;
        end
        %}
        trace_footprint(squeeze(footprint_mtx{isesh}(current_cell,:,:)), plot_color, spc(isesh,1), spc(isesh,2), linewidth);

    end
    
    % probe 8
    if cell_regist_mtx(icell,20)~=0
        current_cell = cell_regist_mtx(icell,20);
        
        if terminal_probe_num>1 && terminal_probe_num<=7 && sum(cell_regist_mtx(icell,terminal_probe_sesh-2:terminal_probe_sesh-1))>0 && active_during_terminal==1
            %plot_color = all_bcolors(terminal_probe_num, :);
            plot_color = all_bcolors(origin_sesh, :);
            trace_footprint(squeeze(footprint_mtx{20}(current_cell,:,:)), plot_color, spc(20,1), spc(20,2), 2);
        elseif ismember(origin_sesh, probe_sessions(2:end-1))
            plot_color = all_bcolors(origin_sesh, :);
            trace_footprint(squeeze(footprint_mtx{20}(current_cell,:,:)), plot_color, spc(20,1), spc(20,2), 1);
        else
            plot_color = .8.*[1 1 1];
            trace_footprint(squeeze(footprint_mtx{20}(current_cell,:,:)), plot_color, spc(20,1), spc(20,2), 1);
        end
        
    end
    
    
end
%}

% aesthetics
axis normal
axis auto
set(gca,'TickLength',[0, 0]); box off;
set(gca, 'YDir','normal')
axis off
set(gcf, 'Position', [217         881        1786         428])
title('Bottom up')


% all origin sessions
origin_sesh_vect = nan(size(cell_regist_mtx,1),1);
for icell = 1:size(cell_regist_mtx,1)
   origin_sesh_vect(icell) = find(cell_regist_mtx(icell,:)>0, 1, 'first'); 
end

% cell regist mtx with only cells that originated during that problem
prob_session = reshape(1:18, 3, 6); prob_session = prob_session(2:3,:);
problem_crms = cell(1,6);
problem_origin_ct = nan(1,6);
for iproblem = 1:6
    problem_crms{iproblem} = cell_regist_mtx(ismember(origin_sesh_vect, prob_session(:, iproblem)),:);
    problem_origin_ct(iproblem) = size(problem_crms{iproblem},1);
end

% compute number of bottom up cells survive to the probe
bottom_up_cells = zeros(6,2);
for iproblem = 1:6
    active_during_next_probe_idx = problem_crms{iproblem}(:,max(prob_session(:, iproblem))+1)>0;
    bottom_up_cells(iproblem, 1) = sum(problem_crms{iproblem}(:,min(prob_session(:, iproblem)))>0 & active_during_next_probe_idx);
    bottom_up_cells(iproblem, 2) = sum(problem_crms{iproblem}(:,max(prob_session(:, iproblem)))>0 & active_during_next_probe_idx);
end

% cells originating in each probe
all_active_probes = nan(6,1);
for iprobe = 1:8
    all_active_probes(iprobe) = sum(cell_regist_mtx(:,probe_sessions(iprobe))>0);
end

% all cells originating in each probe 2:7
probe_originating = nan(6,1);
for iprobe = 2:7
    po_idx = sum(cell_regist_mtx(:,1:probe_sessions(iprobe)-1),2)==0 & cell_regist_mtx(:,probe_sessions(iprobe))>0;
    probe_originating(iprobe-1) = sum( po_idx );
end

% all active cells in each problem
all_active_problems = nan(6,1);
for iprob = 1:6
    all_active_problems(iprob) = sum(sum(cell_regist_mtx(:,prob_session(:,iprob)),2)>0);
end



%% proportion of bottom up cells active in each probe session
figure;hold on

problem_sessions = reshape(problem_sessions,2,size(problem_sessions,2)/2);
problem_sessions = problem_sessions(1,:);
iprob_bar_input = [bottom_up_cells probe_originating]./all_active_probes(2:7);
for iprob = 1:6

    ba = bar([iprob nan], [iprob_bar_input(iprob,:); nan(size(iprob_bar_input(iprob,:)))], 'stacked', 'FaceColor', 'flat');
        session_ct = problem_sessions(iprob);
        for icolor = 1:size(iprob_bar_input,2)
            ba(icolor).CData = all_bcolors(session_ct,:);
            session_ct = session_ct+1;
        end
end
set(gca,'TickLength',[0, 0]); box off;
title('Origin of probe neurons')
xlabel('Probe')
ylabel('Proportion of probe neurons')
%xlim([.25 6.75])
xticks(1:6)
xticklabels(2:7)
ylim([0 1])
set(gcf, 'Position', [1316 851 412 420])


%% proportion of active problem cells that originated in that problem
%{
newcells_problem_prop = problem_origin_ct'./all_active_problems(1:6);
figure;
bar(newcells_problem_prop);
set(gca,'TickLength',[0, 0]); box off;
title('Neurons originating during each problem')
xticklabels(1:6)
xlabel('Problem number')
ylabel('Proportion of active neurons')
ylim([0 1])
%}


%% proportion of problem cells that are reactivated in next probe
%{
bottom_up_cells_problem_prop = bottom_up_cells./all_active_problems(1:6);
figure;
bar(bottom_up_cells_problem_prop);
set(gca,'TickLength',[0, 0]); box off;
title('Active that are reactivated in next probe')
xticklabels(1:6)
xlabel('Problem number')
ylabel('Proportion of active neurons')
ylim([0 1])
%}


%% proportion of new problem cells that are reactivated in next probe
%{
bottom_up_cells_problem_prop = sum(bottom_up_cells,2)./problem_origin_ct';
figure;
bar(bottom_up_cells_problem_prop);
set(gca,'TickLength',[0, 0]); box off;
title('Problem-originating that are reactivated in next probe')
xticklabels(1:6)
xlabel('Problem number')
ylabel('Proportion of new neurons')
ylim([0 1])
set(gcf, 'Position', [1316         851         412         420])
%}


%% proportion of active p8 neurons that are bottom up cells from each probe
% proportion of bottom up cells from each probe that are also active in p8
bottom_up_cells_problem1_p8 = nan(6,1);
bottom_up_cells_problem2_p8 = nan(6,1);
probe_originating_p8 = nan(6,1);
for iproblem = 1:6
    
    active_during_next_probe_idx = cell_regist_mtx(:,max(prob_session(:, iproblem))+1)>0;
    
    % num first problem cells that are active in subsequent probe and p8
    first_problem_origin_idx = origin_sesh_vect==prob_session(1, iproblem);
    bottom_up_cells_problem1_p8(iproblem) = sum(cell_regist_mtx(first_problem_origin_idx & active_during_next_probe_idx, 20)>0);
    
    % num second problem cells that are active in subsequent probe and p8
    second_problem_origin_idx = origin_sesh_vect==prob_session(2, iproblem);
    bottom_up_cells_problem2_p8(iproblem) = sum(cell_regist_mtx(second_problem_origin_idx & active_during_next_probe_idx, 20)>0);
    
    % num subsequent probe originating cells that are also active in p8
    probe_origins_idx = sum(cell_regist_mtx(:,1:probe_sessions(iproblem)),2)==0 & cell_regist_mtx(:,probe_sessions(iproblem+1))>0;
    probe_originating_p8(iproblem) = sum(cell_regist_mtx(probe_origins_idx, 20)>0);
    
end

% bar plot
iprob_bar_input = [bottom_up_cells_problem1_p8 bottom_up_cells_problem2_p8 probe_originating_p8]./sum(cell_regist_mtx(:,probe_sessions(end))>0);
figure; hold on
for iprob = 1:6
    ba = bar([iprob nan], [iprob_bar_input(iprob,:); nan(size(iprob_bar_input(iprob,:)))], 'stacked', 'FaceColor', 'flat');
        session_ct = problem_sessions(iprob);
        for icolor = 1:size(iprob_bar_input,2)
            ba(icolor).CData = all_bcolors(session_ct,:);
            session_ct = session_ct+1;
        end
end


set(gca,'TickLength',[0, 0]); box off;
title('Origin of probe 8 neurons')
xlabel('Probe')
ylabel('Proportion of probe 8 neurons')
ylim([0 1])
%xlim([.25 6.75])
set(gcf, 'Position', [1316         851         412         420])






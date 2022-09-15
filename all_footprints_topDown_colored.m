function [top_down_idx] = all_footprints_topDown_colored(footprint_mtx, cell_regist_mtx)
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
probe_sessions_TDcandidates = probe_sessions(1:6);
session_list = reshape(1:18, 3, 6);

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


%% footprints plot


% top down index
top_down_idx = zeros(size(cell_regist_mtx,1), 6);
for icell = 1:size(cell_regist_mtx,1)
    for iprobe = 1:size(session_list,2)
        % if active on this probe, and subsequent problem, and final probe8
        problem_sessions = probe_sessions(iprobe)+1:probe_sessions(iprobe)+2;
        if cell_regist_mtx(icell, probe_sessions(iprobe))>0 && sum(cell_regist_mtx(icell,problem_sessions),2)>0 && cell_regist_mtx(icell, probe_sessions(end))>0
            top_down_idx(icell, iprobe) = 2;
        elseif cell_regist_mtx(icell, probe_sessions(iprobe))>0 && sum(cell_regist_mtx(icell,problem_sessions),2)>0
            top_down_idx(icell, iprobe) = 1;
        end
    end
end
%top_down_idx_out = zeros(size(cell_regist_mtx));
%top_down_idx_out(:,1:3:13) = top_down_idx;
%
historical_cell_cts = zeros(20, 6);
for icell = 1:size(cell_regist_mtx,1)
    for iprobe = 1:size(session_list,2)
        for isesh = session_list(:,iprobe)'

            current_cell = cell_regist_mtx(icell,isesh);

            if current_cell==0
                continue
            end

            problem_sessions = probe_sessions(iprobe)+1:probe_sessions(iprobe)+2;
            if cell_regist_mtx(icell, probe_sessions(iprobe))>0 && sum(cell_regist_mtx(icell, problem_sessions),2)>0
                % CHANGE COLORING
                %plot_color = all_bcolors(probe_sessions(iprobe), :); % LOCAL TOP DOWN 
                origin_sesh = find(cell_regist_mtx(icell,:)>0, 1, 'first');
                plot_color = all_bcolors(origin_sesh, :); % HISTORICAL TOP DOWN
                historical_cell_cts(origin_sesh, iprobe) = historical_cell_cts(origin_sesh, iprobe)+1;
                
            else
                plot_color = .8.*[1 1 1];
            end
            
            linewidth=1;
            %{
            if cell_regist_mtx(icell, probe_sessions(iprobe))>0 && sum(cell_regist_mtx(icell, problem_sessions),2)>0
                linewidth = 2;
            else
                linewidth = 1;
            end
            %}
            % plot
            trace_footprint(squeeze(footprint_mtx{isesh}(current_cell,:,:)), plot_color, spc(isesh,1), spc(isesh,2), linewidth);

                        
        end
    end
    
    % probe 7
    if cell_regist_mtx(icell,19)~=0
        plot_color = .8.*[1 1 1];
        current_cell = cell_regist_mtx(icell,19);
        trace_footprint(squeeze(footprint_mtx{19}(current_cell,:,:)), plot_color, spc(19,1), spc(19,2));
    end
    % probe 8
    if cell_regist_mtx(icell,20)~=0 %if active p8
        if sum(top_down_idx(icell,:)==2,2)>0
            %plot_color = all_bcolors(probe_sessions_TDcandidates(find(top_down_idx(icell,:)==2 >0, 1, 'first')),:);
            plot_color = all_bcolors(find(cell_regist_mtx(icell,:)>0, 1, 'first'),:);
        else
            plot_color = .8.*[1 1 1];
        end
        current_cell = cell_regist_mtx(icell,20);
        trace_footprint(squeeze(footprint_mtx{20}(current_cell,:,:)), plot_color, spc(20,1), spc(20,2));
    end
    
    
end

% aesthetics
axis normal
axis auto
set(gca,'TickLength',[0, 0]); box off;
set(gca, 'YDir','normal')
axis off
set(gcf, 'Position', [217         881        1786         428])
title('Top down')
%}


%% cell regist mtx with only cells that were in probe
probe_crms = cell(1,7);
bottom_up_idx = cell(1,7);
probe_origin_idx = cell(1,7);
top_down_idx_probecrm = cell(1,7);
for iprobe = 1:7
    
    % crm only with cells active during that probe
    probe_crms{iprobe} = cell_regist_mtx(cell_regist_mtx(:, probe_sessions(iprobe))>0,:);
    
    % idx of cells that originated in previous problem (bottom up)
    if iprobe > 1
        problem_sessions = probe_sessions(iprobe)-2:probe_sessions(iprobe)-1;
        bottom_up_idx{iprobe} = nan(size(probe_crms{iprobe},1),1);
        for ineuron = 1:size(probe_crms{iprobe},1)
            bottom_up_idx{iprobe}(ineuron) = ismember(find(probe_crms{iprobe}(ineuron,:)>0, 1, 'first'), problem_sessions);
        end
    end
    
    % idx of cells that originated in this probe (inference)
    probe_origin_idx{iprobe} = nan(size(probe_crms{iprobe},1),1);
    for ineuron = 1:size(probe_crms{iprobe},1)
        probe_origin_idx{iprobe}(ineuron) = ismember(find(probe_crms{iprobe}(ineuron,:)>0, 1, 'first'), probe_sessions(iprobe));
    end
    
    % idx of cells that are also active in next problem (top down)
    if iprobe < 7
        problem_sessions = probe_sessions(iprobe)+1:probe_sessions(iprobe)+2;
        top_down_idx_probecrm{iprobe} = sum(probe_crms{iprobe}(:, problem_sessions),2)>0;
    end
end


%% proportion of top down cells (in each probe) that are bottom up, inference, or other
%{
figure
for iprobe = 1:6

    topDown_and_bottomUp = sum(sum([top_down_idx_probecrm{iprobe} bottom_up_idx{iprobe}],2)==2);
    topDown_and_probeOrigin = sum(sum([top_down_idx_probecrm{iprobe} probe_origin_idx{iprobe}],2)==2);
    proportions = [topDown_and_bottomUp topDown_and_probeOrigin]./sum(top_down_idx_probecrm{iprobe});
    
    subplot(1,6, iprobe)
    bar(proportions)
    ylim([0 1])
    set(gca,'TickLength',[0, 0]); box off;

end
%}



%% proportion of active cells during probe that will be reativated during next problem
 %{
number_of_active_probe_cells = nan(6,1);
number_also_active_during_next_problem = nan(6,1);
proportions = nan(6,1);
for iprobe = 1:6
    
    problem_sessions = probe_sessions(iprobe)+1:probe_sessions(iprobe)+2;
    
    number_of_active_probe_cells(iprobe) = sum(cell_regist_mtx(:, iprobe)>0);    
    
    number_also_active_during_next_problem(iprobe) = sum( cell_regist_mtx(:, iprobe)>0 & sum(cell_regist_mtx(:, problem_sessions),2)>0 );
        
    proportions(iprobe) = number_also_active_during_next_problem(iprobe)./number_of_active_probe_cells(iprobe);
    
end



figure; bar(proportions)
set(gca,'TickLength',[0, 0]); box off;
title('Probe neurons reactivated during next problem')
xticklabels(1:6)
xlabel('Probe number')
ylabel('Proportion of active neurons')
ylim([0 1])
%}


%% proportion of active cells during problem that were active during prior probe
 
number_of_unq_active_problem_cells = nan(6,1);
number_also_active_during_prior_probe = nan(6,1);
number_prior_problem_origin = nan(6,1);
number_prior_probe_origin = nan(6,1);
number_other_origin = nan(6,1);
proportions = nan(6,1);
for iprobe = 1:6
    
    problem_sessions = probe_sessions(iprobe)+1:probe_sessions(iprobe)+2;
        
    
    number_of_unq_active_problem_cells(iprobe) = sum(sum(cell_regist_mtx(:, problem_sessions),2)>0);
    
    active_during_prior_probe_idx = sum(probe_crms{iprobe}(:, problem_sessions),2)>0;
    number_also_active_during_prior_probe(iprobe) = sum(active_during_prior_probe_idx);

        % find origin session of each cell
        origin_sessions = nan(size(probe_crms{iprobe},1),1);
        for ipcell = 1:size(probe_crms{iprobe}, 1)
            origin_sessions(ipcell) = find(probe_crms{iprobe}(ipcell,:), 1, 'first');
        end
        
        number_prior_problem_origin(iprobe) = sum(ismember(origin_sessions, probe_sessions(iprobe)-2:probe_sessions(iprobe)-1) & active_during_prior_probe_idx); 
        number_prior_probe_origin(iprobe) = sum(ismember(origin_sessions, probe_sessions(iprobe))  & active_during_prior_probe_idx);
        number_other_origin(iprobe) = number_also_active_during_prior_probe(iprobe) - (number_prior_problem_origin(iprobe) + number_prior_probe_origin(iprobe));
    
        
    proportions(iprobe) = number_also_active_during_prior_probe(iprobe)./number_of_unq_active_problem_cells(iprobe);
    
end


%figure; bar(proportions)
figure; bar([number_prior_problem_origin number_prior_probe_origin number_other_origin]./number_of_unq_active_problem_cells, 'stacked')

set(gca,'TickLength',[0, 0]); box off;
title('Problem neurons reactivated from prior probe')
xticklabels(1:6)
xlabel('Problem number')
ylabel('Proportion of active neurons')
ylim([0 1])
legend({'prior problem (bottom-up)', 'prior probe', 'other'})
set(gcf, 'Position', [1316         851         412         420])



%%  proportion of top-down neurons that are reactivated in probe 8
number_active_p8_neurons = sum(cell_regist_mtx(:, probe_sessions(end))>0);

% active p8 neurons originating as top-down neurons in each problem
top_down_idx2 = zeros(size(cell_regist_mtx,1),size(session_list,2));
for iprobe = 1:size(session_list,2)
        
    % active during probe
    active_probe = cell_regist_mtx(:, probe_sessions(iprobe))>0;
    
    % active during post-probe problem
    problem_sessions = probe_sessions(iprobe)+1:probe_sessions(iprobe)+2;
    active_problem = sum(cell_regist_mtx(:, problem_sessions),2)>0;
    
    % active during probe 8
    active_p8 = cell_regist_mtx(:, probe_sessions(end))>0;
    
    % all
    top_down_idx2(active_probe & active_problem & active_p8, iprobe) = 1;
end

%figure; subplot(1,2,1); imagesc(top_down_idx); subplot(1,2,2); imagesc(top_down_idx2)

% counts of prop of probe top-down neurons also active in p8
counts_p8_eachprobe = zeros(20,6);
crm_p8 = cell_regist_mtx(cell_regist_mtx(:,end)>0, :);
%for icell = 1:size(cell_regist_mtx,1)
%    first_activity = find(top_down_idx2(icell,:)>0, 1, 'first');
%    counts_p8(first_activity) = counts_p8(first_activity)+1;
%end
for iprobe = 1:6
    
    % isolate top down neurons (active in probe and at least one of the
    % next two (problem) sessions
    probe_sesh = probe_sessions(iprobe);
    top_down_iso = crm_p8(crm_p8(:,probe_sesh)>0,:);
    top_down_iso = top_down_iso(sum(top_down_iso(:,probe_sesh+1:probe_sesh+2),2)>0, :);

    % id origins
    top_down_iso_origins = nan(size(top_down_iso,1),1);
    for icell = 1:size(top_down_iso,1)
        top_down_iso_origins(icell) = find(top_down_iso(icell,:)>0, 1, 'first');
    end
        
    % load 20,1 vector counting them by their origin sessions
    counts_p8_eachprobe(:, iprobe) = histcounts(top_down_iso_origins, 1:20+1)';
    %sum(top_down_idx2(:,iprobe)>0);
end

%counts_p8
bar_input = counts_p8_eachprobe./number_active_p8_neurons;
figure; ba = bar(bar_input', 'stacked', 'FaceColor', 'flat');
for icolor = 1:size(cell_regist_mtx,2)
   ba(icolor).CData = all_bcolors(icolor,:);
end
set(gca,'TickLength',[0, 0]); box off;
title('Origin of probe 8 neurons'); %Proportion of p8 neurons that were top-down neurons in each probe'
xticklabels(1:6)
xlabel('Probe')
ylabel('Proportion of active p8 neurons')
ylim([0 1])
set(gcf, 'Position', [1316 851 412 420])



%% proportion of each problem's top-down neurons that arose in previous probe as either bottom up, new to that probe, or other











function all_footprints_inference_colored(footprint_mtx, cell_regist_mtx)
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

% session plotting coordinates
yspace = 275;
spc = ...
[0 0; 200 -1*yspace; 200 -2*yspace; ...
400 0; 600 -1*yspace; 600 -2*yspace; ...
800 0; 1000 -1*yspace; 1000 -2*yspace; ...
1200 0; 1400 -1*yspace; 1400 -2*yspace; ...
1600 0; 1800 -1*yspace; 1800 -2*yspace; ...
2000 0; 2200 -1*yspace; 2200 -2*yspace; ...
2400 0; 2800 0];


% plot
%
for icell = 1:size(cell_regist_mtx,1)
    
    % first session for this cell (term sessions will only plot if cell
    % active)
    origin_sesh = find(cell_regist_mtx(icell,:)>0, 1, 'first');

    for isesh = 1:size(cell_regist_mtx,2)-1

        current_cell = cell_regist_mtx(icell,isesh);

        if current_cell==0
            continue
        end

        if ismember(origin_sesh, probe_sessions) && origin_sesh==isesh
            plot_color = all_bcolors(probe_sessions(ismember(probe_sessions, origin_sesh)), :);
        else
            plot_color = .8.*[1 1 1];
        end
        
        linewidth = 1;
        
        %isesh
        %current_cell
        %spc(isesh,1);
        %spc(isesh,2);
        %squeeze(footprint_mtx{isesh}(current_cell,:,:));

        trace_footprint(squeeze(footprint_mtx{isesh}(current_cell,:,:)), plot_color, spc(isesh,1), spc(isesh,2), linewidth);

    end
    
    % probe 8
    if cell_regist_mtx(icell,20)>0
        current_cell = cell_regist_mtx(icell,20);
        if ismember(origin_sesh, probe_sessions(1:7))
            plot_color = all_bcolors(probe_sessions(ismember(probe_sessions, find(cell_regist_mtx(icell,:)>0, 1, 'first'))), :);
        else
            plot_color = .8.*[1 1 1];
        end
        trace_footprint(squeeze(footprint_mtx{20}(current_cell,:,:)), plot_color, spc(20,1), spc(20,2));
    end

end
%}

% aesthetics
axis normal
axis on
axis auto
set(gca,'TickLength',[0, 0]); box off;
set(gca, 'YDir','normal')
set(gcf, 'Position', [217 881 1786 428])
title('Probe-originating')


% all origin sessions
origin_sesh_vect = nan(size(cell_regist_mtx,1),1);
for icell = 1:size(cell_regist_mtx,1)
   origin_sesh_vect(icell) = find(cell_regist_mtx(icell,:)>0, 1, 'first'); 
end

% number of cells originating during each probe
probe_origin_cts = histcounts(origin_sesh_vect, 1:size(cell_regist_mtx)+realmin);
probe_origin_cts = probe_origin_cts(probe_sessions);

% cell regist mtx with only cells that originated during that probe
probe_crms = cell(1,7);
for iprobe = 1:7
    probe_crms{iprobe} = cell_regist_mtx(ismember(origin_sesh_vect, probe_sessions(iprobe)),:);
end

% all active cells in each probe 1:7
all_active_probes = nan(6,1);
for iprobe = 1:8
    all_active_probes(iprobe) = sum(cell_regist_mtx(:,probe_sessions(iprobe))>0);
end






%% proportion of active problem cells that originated in that problem
newcells_probe_prop = probe_origin_cts(1:7)'./all_active_probes(1:7);
figure;
bar(newcells_probe_prop);
set(gca,'TickLength',[0, 0]); box off;
title('Neurons originating during each probe')
xticklabels(1:6)
xlabel('Probe number')
ylabel('Proportion of active neurons')
ylim([0 1])



%% compute proportion p8 cells that originated from each probe
active_during_p8_ct = zeros(6,1);
for iprobe = 1:6
    active_during_p8_ct(iprobe) = sum(probe_crms{iprobe}(:,probe_sessions(end))>0);
end
figure;
bar(active_during_p8_ct./all_active_probes(8))
set(gca,'TickLength',[0, 0]); box off;
title('Probe-originating neurons reactivated in probe 8')
xlabel('Originating problem')
ylabel('Proportion of active p8 neurons')
ylim([0 1])








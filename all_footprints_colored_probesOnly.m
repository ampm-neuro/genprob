function all_footprints_colored_probesOnly(footprint_mtx, cell_regist_mtx)
% plots footprints for all probes colored according to which probe
%
% try aligned_data_struct.spatial_footprints_corrected for spatial
% footprints
%

figure; hold on

% colors
%
all_bcolors = [.3.*[255 255 255];...
    186 41 41; ...
    232 116 0; ...
    232 204 37; ...
    28 142 22; ...
    22 178 247; ...
    48 83 193; ....
    143 020 181]./255;
%}

% session order
session_order = [1 9 15 2 10 16 3 11 17 4 12 18 5 13 19 6 14 20 7 8];
cell_regist_mtx = cell_regist_mtx(:,session_order);
footprint_mtx = footprint_mtx(session_order);

% probe sessions only
probe_sessions = [1:3:19 20];
cell_regist_mtx = cell_regist_mtx(:, probe_sessions);
footprint_mtx = footprint_mtx(probe_sessions);


% session plotting coordinates
spc_incr = 275;
spc = 0:spc_incr:spc_incr*size(cell_regist_mtx,2);


% plot
for icell = 1:size(cell_regist_mtx,1)
        
    plot_color = all_bcolors(find(cell_regist_mtx(icell,:)>0, 1, 'first'),:);
    
    for isesh = 1:size(cell_regist_mtx,2)
        
        current_cell = cell_regist_mtx(icell,isesh);
        
        if current_cell==0
            continue
        end
             
        trace_footprint(squeeze(footprint_mtx{isesh}(current_cell,:,:)), plot_color, spc(isesh));

    end
    
end


% aesthetics
axis normal
axis on
axis auto
set(gca,'TickLength',[0, 0]); box off;
set(gca, 'YDir','normal')
set(gcf, 'Position', [217         881        1786         428])
title('Unique session color')



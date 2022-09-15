function all_footprints_quads_colored(footprint_mtx, cell_regist_mtx)
% plots footprints for all probes colored according to which probe
%
% try aligned_data_struct.spatial_footprints_corrected for spatial
% footprints
%

% colors
all_bcolors = [097 191 050; ...
               186 041 041; ...
               051 153 255; ...
               255 128 000]; ...
               %143 020 181; ...
               %249 221 037; ...
               %204 000 102; ...
               %120 120 120];
all_bcolors = all_bcolors./255;
%}


% session order
session_order = [1 9 15 2 10 16 3 11 17 4 12 18 5 13 19 6 14 20 7 8];
cell_regist_mtx = cell_regist_mtx(:,session_order);
footprint_mtx = footprint_mtx(session_order);

% quads
quads = [1 2 3 4; 4 5 6 7; 7 8 9 10; 10 11 12 13; 13 14 15 16; 16 17 18 19];

subplot_idx = [1:4; 9:12; 17:20; 5:8; 13:16; 21:24];

figure; hold on

% iterate through quads
%{
for iquad = 1:size(quads,1)
    % local vars
    local_crm = cell_regist_mtx(:,quads(iquad,:));
    local_fm = footprint_mtx(quads(iquad,:));
    
    % plot
    subplot(3, 8, subplot_idx(iquad,:))
    for ineuron = 1:size(local_crm,1)

        plot_color = all_bcolors(find(local_crm(ineuron,:)>0, 1, 'first'),:);

        isesh_ct = 0;
        mean_xpos = nan(1,size(local_crm,2));
        for isesh = 1:size(local_crm,2)
            isesh_ct = isesh_ct +1;
            current_cell = local_crm(ineuron,isesh);

            if current_cell==0
                continue
            end
            
            mean_xpos(isesh) = isesh_ct*300;
            trace_footprint(squeeze(local_fm{isesh}(current_cell,:,:)), plot_color, mean_xpos(isesh));

        end

    end

    % aesthetics
    axis normal
    axis on
    axis auto
    set(gca,'TickLength',[0, 0]); box off;
    title(['Problem 0' num2str(iquad)])
    xticklabels({'Probe Before', 'Problem first', 'Problem last', 'Probe After'})
    
end
%}

figure; hold on

% bar plot of proportions
for iquad = 1:size(quads,1)
    
    % local vars
    local_crm = cell_regist_mtx(:,quads(iquad,:));
    origin_mtx = nan(size(local_crm));
    
    for ineuron = 1:size(local_crm,1)
        origin_mtx(ineuron, local_crm(ineuron,:)>0) = find(local_crm(ineuron,:)>0, 1, 'first');
    end
    
    % proportion plots
    subplot(3, 8, subplot_idx(iquad,:)); hold on
    mean_xpos = nan(1,size(origin_mtx,2));
    proportions = nan(size(origin_mtx,2), 4);
    for isesh = 1:size(origin_mtx,2)
        
        first_sesh_cells = sum(origin_mtx(:,isesh)==1);
        second_sesh_cells = sum(origin_mtx(:,isesh)==2);
        third_sesh_cells = sum(origin_mtx(:,isesh)==3);
        fourth_sesh_cells = sum(origin_mtx(:,isesh)==4);
        proportions(isesh,:) = [first_sesh_cells second_sesh_cells third_sesh_cells fourth_sesh_cells]./sum(~isnan(origin_mtx(:,isesh)));
        
        xpos = quads(isesh,:)+(isesh-1);
        mean_xpos(isesh) = mean(xpos);
        
    end
    
    hb = bar(mean_xpos, proportions);
    set(gca,'TickLength',[0, 0]); box off;
    
    hb(1).FaceColor = all_bcolors(1,:);
    hb(2).FaceColor = all_bcolors(2,:);
    hb(3).FaceColor = all_bcolors(3,:);
    hb(4).FaceColor = all_bcolors(4,:);
    
    xticks(mean_xpos)
    xticklabels({'Probe Before', 'Problem first', 'Problem last', 'Probe After'})
    
    if ismember(iquad ,1:3)
        ylabel('Proportion')
        yticks([0 1])
    else
        yticks([])
    end
end




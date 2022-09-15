green_color = [21 154 72]./255; blue_color =  [43 73 159]./255;
light_green_color = [162 201 129]./255; light_blue_color = [102 164 213]./255;
figure_path = 'C:\Users\ampm1\Documents\manuscripts\gentrain\manuscript_figures\sfig03\';
predict_subjs = [{'651049m1'} {'658648m2'} {'691359m1'} {'690330m1'} {'691359m2'} {'696944m3'} {'699437m4'} {'699438m3'}];
unpredict_subjs = [{'687034m4'} {'683472m2'} {'699437m2'} {'699438m1'} {'683472m3'}];

% for time warping
tses = [0.2 1.0 1.0 2.0 3.0 2.0];


%% Fields of view and footprints rows
%{
folder_path = 'G:\683472m2_tozip';
plot_overlay_footprints_probes_overlay(folder_path)
var_name = 'FOV_footprints_rows'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
%}

for isesh = 1:7
    load(['C:\Users\ampm1\Desktop\sfig3temp\sfig3temp\p' num2str(isesh) '\intermediate_results.mat'], 'initialization'); 
    figure; imagesc(initialization.neuron.Cn); colormap(bone)
    title(['Probe ' num2str(isesh-1)])
    axis square; axis off
    var_name = ['MaxProj_' num2str(isesh)]; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
end

%% adjacent probe overlapping cells

load('C:\Users\ampm1\Documents\MATLAB\generalization\cellreg_data\696944m3\cell_reg_696944m3.mat');
probe_sessions = 1:7;
for isesh = probe_sessions % session 1 (and 2) was used in figure
    figure; hold on
    trace_footprints_sessions_overlap(cell_registered_struct.spatial_footprints_corrected([isesh isesh+1]), cell_registered_struct.cell_to_index_map(:,[isesh isesh+1]))
    legend('location', 'northeastoutside')
    legend({'session 2', 'session 1'})
    close; % close venn diagram
    var_name = ['Overlapping_population_example_footprints_' num2str(isesh)]; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
end


%% example reliable firing fields in held cells
% plot all trials (across sessions) from single neurons

% only chosen cell numbers.. or comment out
chosen_cell_nums = [2 17];

% iterate through subjects
min_active_sessions = 7;
for imouse = {'696944m3'} %[predict_subjs unpredict_subjs] 

    mouse = imouse{1}

load(['C:\Users\ampm1\Documents\MATLAB\generalization\cellreg_data\' mouse '\cell_reg_' mouse '.mat']); 

% cell registration
cell_regist_mtx = cell_registered_struct.cell_to_index_map; % cell_regist_mtx
cell_regist_mtx_sorted = cell_regist_mtx(:,1:7); % probes only

    % filter
    max_num_active_sessions = max(sum(cell_regist_mtx_sorted>0,2));
    if max_num_active_sessions < min_active_sessions
        continue
    end
    
% spatial footprints
footprint_mtx = cell_registered_struct.spatial_footprints_corrected;
footprint_mtx = footprint_mtx(1:7); % probes only

% probe sessions 0-6
session_cell = [...
get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_imaging_hpc\' mouse ], 'preprobe', 'LED');...
get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_imaging_hpc\' mouse ] , 'postprobe', 'LED');...
get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_hivar_imaging_hpc\' mouse ], 'preprobe', 'LED');...
get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_hivar_imaging_hpc\' mouse ] , 'postprobe', 'LED')];
session_cell_sorted = session_cell(1:7);

% activity from every cell
[~, session_mtx_cell] = image_subject_allsesh_decode(session_cell_sorted, cell_regist_mtx_sorted, tses);

% ony cells active on all probe sessions
cells_with_min_active_sessions = find(sum(cell_regist_mtx_sorted>0,2)>=min_active_sessions);
for ic = cells_with_min_active_sessions'
    
    % only chosen egs
    if exist(chosen_cell_nums, 'var') && ~ismember(ic, [chosen_cell_nums])
        continue
    end
    
    all_trials_one_cell_fr(ic, session_mtx_cell, tses);
    title([mouse '; cell ' num2str(ic)])
    set(gcf, 'Position', [82 377 600 500])
    var_name = 'image_allprobes_eg_1'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
    
    % add footprint tracking
    %
    figure;
    for isesh = 1:size(session_cell_sorted,1) 
        disp(num2str(cell_regist_mtx_sorted(ic, isesh))); 
        subplot(1,size(session_cell_sorted,1),isesh); 
        trace_footprints(footprint_mtx{isesh}, 0.8.*[1 1 1]); 
        hold on; 
        trace_footprint(squeeze(footprint_mtx{isesh}(cell_regist_mtx_sorted(ic, isesh),:,:)), [0 0 0]);
    end
    sgtitle([mouse '; cell ' num2str(ic)])
    set(gcf, 'Position', [704 380 1088 500])
    var_name = 'image_allprobes_eg_1_footprints'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')
    %}
    
end
    
    

end

%% Spatial correlations vs centroid (anatomical) distances

% all cellreg folders
subject_folders = get_folder_paths_all('C:\Users\ampm1\Documents\MATLAB\generalization\cellreg_data', 0);

% iterate through predict subjects
%
for isubj = 1%:length(predict_subjs)
   isubj_predict = isubj
    % subject data
    spatial_correlations = [];
    centroid_distances = [];
    info_output_mtx = [];
    overall_cell = 0;
    
    % subject folder
    subj_cr_folder = subject_folders{contains(subject_folders, predict_subjs{isubj})};
    
    % load cell registration structure
    cell_reg_path = get_file_paths_targeted(subj_cr_folder, {'cell_reg', '.mat'});
    load(cell_reg_path{1}, 'cell_registered_struct');
    
    % iterate through probe session pairs
    for isesh1 = 1:8
        for isesh2 = 1:8
            tic
            %[isesh1 isesh2]
            
            % skip duplicates and identity
            if isesh2 <= isesh1; continue; end

            % spatial footprint data
            sp_sesh1 = cell_registered_struct.spatial_footprints_corrected{isesh1};
            sp_sesh2 = cell_registered_struct.spatial_footprints_corrected{isesh2};
            % centroid position data
            cp_sesh1 = cell_registered_struct.centroid_locations_corrected{isesh1};
            cp_sesh2 = cell_registered_struct.centroid_locations_corrected{isesh2};
            % registration
            cell_reg_mtx = cell_registered_struct.cell_to_index_map(:,[isesh1 isesh2]);

            % add rows for every combination of neurons
            neuron_combos = size(sp_sesh1,1)*size(sp_sesh2,1)
            spatial_correlations = [spatial_correlations; nan(neuron_combos, 1)];
            centroid_distances = [centroid_distances; nan(neuron_combos, 1)];
            info_output_mtx = [info_output_mtx; nan(neuron_combos, 6)];

            % iterate through neurons pairs
            for ic1 = 1:size(sp_sesh1,1)
                for ic2 = 1:size(sp_sesh2,1)
                    overall_cell = overall_cell + 1;

                    % compute pairwise spatial footprint correlation
                    %
                        % cell size and position
                        %halfedge = 15*1.3; % 1.3 microns per pixel
                        %round_centroid_c1_s1 = round(cp_sesh1(ic1,:))
                        %round_centroid_c2_s2 = round(cp_sesh2(ic2,:))

                        % indices
                        %xy_c1_s1 = round([cp_sesh1(ic1,2)-halfedge, cp_sesh1(ic1,2)+halfedge, cp_sesh1(ic1,1)-halfedge, cp_sesh1(ic1,1)+halfedge]);
                        %xy_c2_s2 = round([cp_sesh2(ic2,2)-halfedge, cp_sesh2(ic2,2)+halfedge, cp_sesh2(ic2,1)-halfedge, cp_sesh2(ic2,1)+halfedge]);

                        % skip cell if mask indices extend beyond FoV
                        %{
                        mask_c1_sesh1 = squeeze(sp_sesh1(ic1,:,:));
                        mask_c2_sesh2 = squeeze(sp_sesh2(ic2,:,:));
                        if xy_c1_s1(1)<1 || xy_c1_s1(2)>size(mask_c1_sesh1,1) || xy_c1_s1(3)<1 || xy_c1_s1(4)>size(mask_c1_sesh1,2)
                            continue
                        elseif xy_c2_s2(1)<1 || xy_c2_s2(2)>size(mask_c2_sesh2,1) || xy_c2_s2(3)<1 || xy_c2_s2(4)>size(mask_c2_sesh2,2)
                            continue
                        end
%}
                        % index for masks
                        %mask_c1_sesh1 = mask_c1_sesh1(xy_c1_s1(1):xy_c1_s1(2), xy_c1_s1(3):xy_c1_s1(4));                        
                        %mask_c2_sesh2 = mask_c2_sesh2(xy_c2_s2(1):xy_c2_s2(2), xy_c2_s2(3):xy_c2_s2(4));
                    
                        % 
                                                
                        % binarize
                        %mask_c1_sesh1(mask_c1_sesh1>0) = 1;
                        %mask_c2_sesh2(mask_c2_sesh2>0) = 1;
                        
                        % vectorize
                        c1_sesh1 = sp_sesh1(ic1,:,:);
                        c2_sesh2 = sp_sesh2(ic2,:,:);
                        
                        % spatial correlation
                        spatial_corr = corr(c1_sesh1(:), c2_sesh2(:));
                        %spatial_corr = corr(mask_c1_sesh1(:), mask_c2_sesh2(:));
                        
                        
                        % load spatial correlation
                        spatial_correlations(overall_cell) = spatial_corr;

                    % compute anatomical centroid distance
                    %
                        % distance
                        centroid_distance = pdist([cp_sesh1(ic1,:); cp_sesh2(ic2,:)]);
                        % load distance
                        centroid_distances(overall_cell) = centroid_distance;

                    % registration
                    %    
                        % registered as same cell test
                        registered_SameOrDif = cell_reg_mtx(cell_reg_mtx(:,1)==ic1, 2)==ic2; % 1 for same cell, 0 for different cells
                        if isempty(registered_SameOrDif)
                            registered_SameOrDif = nan;
                        end

                        % load
                        info_output_mtx(overall_cell,:) =  [isubj isesh1 ic1 isesh2 ic2 registered_SameOrDif];
                end
            end
            toc/neuron_combos
        end
    end
    
    % save subject data
    save(['SFigure03_predict_' num2str(isubj)], 'centroid_distances', 'spatial_correlations', 'info_output_mtx');    
    
end
%}



% iterate through unpredict subjects
%
for isubj = 1:length(unpredict_subjs)-1
   isubj_unpredict = isubj
    % subject data
    spatial_correlations = [];
    centroid_distances = [];
    info_output_mtx = [];
    overall_cell = 0;
    
    % subject folder
    subj_cr_folder = subject_folders{contains(subject_folders, unpredict_subjs{isubj})};
    
    % load cell registration structure
    cell_reg_path = get_file_paths_targeted(subj_cr_folder, {'cell_reg', '.mat'});
    load(cell_reg_path{1}, 'cell_registered_struct');
    
    % iterate through probe session pairs
    for isesh1 = 1:8
        for isesh2 = 1:8
            tic
            %[isesh1]
            
            % skip duplicates and identity
            if isesh2 <= isesh1; continue; end

            % spatial footprint data
            sp_sesh1 = cell_registered_struct.spatial_footprints_corrected{isesh1};
            sp_sesh2 = cell_registered_struct.spatial_footprints_corrected{isesh2};
            % centroid position data
            cp_sesh1 = cell_registered_struct.centroid_locations_corrected{isesh1};
            cp_sesh2 = cell_registered_struct.centroid_locations_corrected{isesh2};
            % registration
            cell_reg_mtx = cell_registered_struct.cell_to_index_map(:,[isesh1 isesh2]);

            % add rows for every combination of neurons
            neuron_combos = size(sp_sesh1,1)*size(sp_sesh2,1)
            spatial_correlations = [spatial_correlations; nan(neuron_combos, 1)];
            centroid_distances = [centroid_distances; nan(neuron_combos, 1)];
            info_output_mtx = [info_output_mtx; nan(neuron_combos, 6)];

            % iterate through neurons pairs
            for ic1 = 1:size(sp_sesh1,1)
                for ic2 = 1:size(sp_sesh2,1)
                    overall_cell = overall_cell + 1;

                    % compute pairwise spatial footprint correlation
                    %
                        % cell size and position
                        %{
                        halfedge = 15*1.3; % 1.3 microns per pixel
                        round_centroid_c1_s1 = round(cp_sesh1(ic1,:));
                        round_centroid_c2_s2 = round(cp_sesh2(ic2,:));

                        % indices
                        xy_c1_s1 = [round_centroid_c1_s1(2)-halfedge, round_centroid_c1_s1(2)+halfedge, round_centroid_c1_s1(1)-halfedge, round_centroid_c1_s1(1)+halfedge];
                        xy_c2_s2 = [round_centroid_c2_s2(2)-halfedge, round_centroid_c2_s2(2)+halfedge, round_centroid_c2_s2(1)-halfedge, round_centroid_c2_s2(1)+halfedge];

                        % skip cell if mask indices extend beyond FoV
                        mask_c1_sesh1 = squeeze(sp_sesh1(ic1,:,:));
                        mask_c2_sesh2 = squeeze(sp_sesh2(ic2,:,:));
                        if xy_c1_s1(1)<1 || xy_c1_s1(2)>size(mask_c1_sesh1,1) || xy_c1_s1(3)<1 || xy_c1_s1(4)>size(mask_c1_sesh1,2)
                            continue
                        elseif xy_c2_s2(1)<1 || xy_c2_s2(2)>size(mask_c2_sesh2,1) || xy_c2_s2(3)<1 || xy_c2_s2(4)>size(mask_c2_sesh2,2)
                            continue
                        end

                        % index for masks
                        mask_c1_sesh1 = mask_c1_sesh1(xy_c1_s1(1):xy_c1_s1(2), xy_c1_s1(3):xy_c1_s1(4));                        
                        mask_c2_sesh2 = mask_c2_sesh2(xy_c2_s2(1):xy_c2_s2(2), xy_c2_s2(3):xy_c2_s2(4));
                        
                        % spatial correlation
                        spatial_corr = corr(mask_c1_sesh1(:), mask_c2_sesh2(:));
                        %}
                        % vectorize
                        c1_sesh1 = sp_sesh1(ic1,:,:);
                        c2_sesh2 = sp_sesh2(ic2,:,:);
                        
                        % spatial correlation
                        spatial_corr = corr(c1_sesh1(:), c2_sesh2(:));
                        
                    % compute anatomical centroid distance
                    %
                        % distance
                        centroid_distance = pdist([cp_sesh1(ic1,:); cp_sesh2(ic2,:)]);
                        % load distance
                        centroid_distances(overall_cell) = centroid_distance;

                    % registration
                    %    
                        % registered as same cell test
                        registered_SameOrDif = cell_reg_mtx(cell_reg_mtx(:,1)==ic1, 2)==ic2; % 1 for same cell, 0 for different cells
                        if isempty(registered_SameOrDif)
                            registered_SameOrDif = nan;
                        end

                        % load
                        info_output_mtx(overall_cell,:) =  [isubj isesh1 ic1 isesh2 ic2 registered_SameOrDif];
                end
            end
            toc/neuron_combos
        end
    end
    
    % save subject data
    save(['SFigure03_unpredict_' num2str(isubj)], 'centroid_distances', 'spatial_correlations', 'info_output_mtx');    
    
end
%}

%% merge predictable subjects
spatial_correlations_predict = [];
centroid_distances_predict = [];
info_output_mtx_predict = [];
for isubj = 1%:length(predict_subjs)
    load(['SFigure03_predict_' num2str(isubj) '.mat'], 'centroid_distances', 'spatial_correlations', 'info_output_mtx');
    spatial_correlations_predict = [spatial_correlations_predict; spatial_correlations];
    centroid_distances_predict = [centroid_distances_predict; centroid_distances];
    info_output_mtx_predict = [info_output_mtx_predict; info_output_mtx];
end

% merge unpredictable subjects
spatial_correlations_unpredict = [];
centroid_distances_unpredict = [];
info_output_mtx_unpredict = [];
for isubj = 1:length(unpredict_subjs)-1
    load(['SFigure03_unpredict_' num2str(isubj) '.mat'], 'centroid_distances', 'spatial_correlations', 'info_output_mtx');
    spatial_correlations_unpredict = [spatial_correlations_unpredict; spatial_correlations];
    centroid_distances_unpredict = [centroid_distances_unpredict; centroid_distances];
    info_output_mtx_unpredict = [info_output_mtx_unpredict; info_output_mtx];
end


%% plot probability heatmap for all neurons in each group
figure; 

subplot(2,2,1)
hist2d(centroid_distances_predict(info_output_mtx_predict(:,6)==1), spatial_correlations_predict(info_output_mtx_predict(:,6)==1), 0:1:26, 0:.05:1, 'tile', 'probability')
set(gca, 'YDir','normal'); axis square; yticks(0:.25:1); ylabel('Spatial correlation (r)'); xticks(0:5.2:26); xticklabels(0:4:20); xlabel('Centroid distance (micron)')
caxis([0 .11]); colorbar
title('Predict, same cell')

subplot(2,2,3)
hist2d(centroid_distances_predict(info_output_mtx_predict(:,6)==0), spatial_correlations_predict(info_output_mtx_predict(:,6)==0), 0:1:26, 0:.05:1, 'tile', 'probability')
set(gca, 'YDir','normal'); axis square; yticks(0:.25:1); ylabel('Spatial correlation (r)'); xticks(0:5.2:26); xticklabels(0:4:20); xlabel('Centroid distance (micron)')
caxis([0 .11]); colorbar
title('Predict, different cells')

subplot(2,2,2)
hist2d(centroid_distances_unpredict(info_output_mtx_unpredict(:,6)==1), spatial_correlations_unpredict(info_output_mtx_unpredict(:,6)==1), 0:1:26, 0:.05:1, 'tile', 'probability')
set(gca, 'YDir','normal'); axis square; yticks(0:.25:1); ylabel('Spatial correlation (r)'); xticks(0:5.2:26); xticklabels(0:4:20); xlabel('Centroid distance (micron)')
caxis([0 .11]); colorbar
title('Unpredict, same cell')

subplot(2,2,4)
hist2d(centroid_distances_unpredict(info_output_mtx_unpredict(:,6)==0), spatial_correlations_unpredict(info_output_mtx_unpredict(:,6)==0), 0:1:26, 0:.05:1, 'tile', 'probability')
set(gca, 'YDir','normal'); axis square; yticks(0:.25:1); ylabel('Spatial correlation (r)'); xticks(0:5.2:26); xticklabels(0:4:20); xlabel('Centroid distance (micron)')
caxis([0 .11]); colorbar
title('Unpredict, different cells')

var_name = 'Heatmap2d_spacecorr_dist'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')

%% plot median spatial correlations for every pairwise session

% predictable
pred_median_mtx_spatial = nan(8,8);
pred_median_mtx_centroid = nan(8,8);
for isesh1 = 1:8
    for isesh2 = 1:8
        if isesh2>isesh1
            pred_median_mtx_spatial(isesh1, isesh2) = median(spatial_correlations_predict(info_output_mtx_predict(:,6)==1 & info_output_mtx_predict(:,2)==isesh1 & info_output_mtx_predict(:,4)==isesh2));
            pred_median_mtx_centroid(isesh1, isesh2) = median(centroid_distances_predict(info_output_mtx_predict(:,6)==1 & info_output_mtx_predict(:,2)==isesh1 & info_output_mtx_predict(:,4)==isesh2));
        else
            pred_median_mtx_spatial(isesh1, isesh2) = median(spatial_correlations_predict(info_output_mtx_predict(:,6)==1 & info_output_mtx_predict(:,2)==isesh2 & info_output_mtx_predict(:,4)==isesh1));
            pred_median_mtx_centroid(isesh1, isesh2) = median(centroid_distances_predict(info_output_mtx_predict(:,6)==1 & info_output_mtx_predict(:,2)==isesh2 & info_output_mtx_predict(:,4)==isesh1));
        end
    end
end

% unpredictable
unpred_median_mtx_spatial = nan(8,8);
unpred_median_mtx_centroid = nan(8,8);
for isesh1 = 1:8
    for isesh2 = 1:8
        if isesh2>isesh1
            unpred_median_mtx_spatial(isesh1, isesh2) = median(spatial_correlations_unpredict(info_output_mtx_unpredict(:,6)==1 & info_output_mtx_unpredict(:,2)==isesh1 & info_output_mtx_unpredict(:,4)==isesh2));
            unpred_median_mtx_centroid(isesh1, isesh2) = median(centroid_distances_unpredict(info_output_mtx_unpredict(:,6)==1 & info_output_mtx_unpredict(:,2)==isesh1 & info_output_mtx_unpredict(:,4)==isesh2));
        else
            unpred_median_mtx_spatial(isesh1, isesh2) = median(spatial_correlations_unpredict(info_output_mtx_unpredict(:,6)==1 & info_output_mtx_unpredict(:,2)==isesh2 & info_output_mtx_unpredict(:,4)==isesh1));
            unpred_median_mtx_centroid(isesh1, isesh2) = median(centroid_distances_unpredict(info_output_mtx_unpredict(:,6)==1 & info_output_mtx_unpredict(:,2)==isesh2 & info_output_mtx_unpredict(:,4)==isesh1));
        end
    end
end

figure
subplot(2,2,1); imagesc(pred_median_mtx_spatial); title('Predict, spatial'); xticks(1:8); yticks(1:8);
axis square; caxis([0.5 1]); colorbar; xlabel('Probe session'); ylabel('Probe session')
subplot(2,2,3); imagesc(pred_median_mtx_centroid); title('Predict, centroid'); xticks(1:8); yticks(1:8);
axis square; caxis([0 1.3*8]); colorbar; xlabel('Probe session'); ylabel('Probe session')
subplot(2,2,2); imagesc(unpred_median_mtx_spatial); title('Unpredict, spatial'); xticks(1:8); yticks(1:8);
axis square; caxis([0.5 1]); colorbar; xlabel('Probe session'); ylabel('Probe session')
subplot(2,2,4); imagesc(unpred_median_mtx_centroid); title('Unpredict, centroid'); xticks(1:8); yticks(1:8);
axis square; caxis([0 1.3*8]); colorbar; xlabel('Probe session'); ylabel('Probe session')

var_name = 'PairwiseProbe_mtx_spacecorr_dist'; print([figure_path var_name], '-dpdf', '-painters', '-bestfit')












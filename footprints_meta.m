function footprints_meta(cell_registered_struct, save_fp)
%combine one footprint per unique neuron into a cell-reg ready footprint
%file
%cell_registered_struct.cell_to_index_map
all_fp_idx = zeros(size(cell_registered_struct.cell_to_index_map,1),1);
footprints_meta_out = nan(size(cell_registered_struct.cell_to_index_map,1),...
    size(cell_registered_struct.spatial_footprints_corrected{1},2),...
    size(cell_registered_struct.spatial_footprints_corrected{1},3));
for isesh = 1:length(cell_registered_struct.spatial_footprints_corrected)
    % in session, not yet in fpm
    sesh_idx = cell_registered_struct.cell_to_index_map(:,isesh)>0 & all_fp_idx==0;
    footprints_meta_out(sesh_idx,:,:) = ...
        cell_registered_struct.spatial_footprints_corrected{isesh}(cell_registered_struct.cell_to_index_map(sesh_idx, isesh),:,:);
end

%save([save_fp '/footprints_meta_out.mat'], 'footprints_meta_out')
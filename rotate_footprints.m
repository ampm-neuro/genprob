function fps = rotate_footprints(fps, cw_deg, translation)
% rotate footprints_in (fps_in) clockwise by ccw_deg degrees

    % hard coded
    microns_per_pixel = .5;
    
    % center of all points
    centroids = nan(size(fps,1),2);
    for ic = 1:size(fps,1)
        local_fp = squeeze(fps(ic,:,:));
        [centroids(ic,1), centroids(ic,2)] = ind2sub(size(local_fp),  find(local_fp==max(local_fp, [], 'all'),1));
    end
    center_of_fov = round(mean(centroids));
    
    % iterate through fps
    for ic = 1:size(fps,1)
        local_fp = squeeze(fps(ic,:,:));

        
        
        %figure; 
        %subplot(1,2,1); imagesc(squeeze(fps(ic,:,:))); title before; caxis_hold = caxis; axis square;
        
        fps(ic,:,:) = rotate_spatial_footprint(local_fp, cw_deg, translation, center_of_fov, centroids(ic,:), microns_per_pixel);
        
        %subplot(1,2,2); imagesc(squeeze(fps(ic,:,:))); title after; caxis(caxis_hold); axis square;
        
        %drawnow
    end
end
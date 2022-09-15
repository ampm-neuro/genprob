function plot_overlay_footprints_probes(folder_path)
% input the folder path containing cnmfe outputs and it will subplot the
% overlay_footprint pdfs from every day, in chronological order

% 'E:\ampm\690330m1'

% sub folders
subfolders = [{'probe01'} {'probe02'} {'probe03'} {'probe04'} {'probe06'} {'probe07'} {'probe08'}];

% figure
figure;
          
% iterate through the sub folders          
iplot = 0;
for isf = 1:length(subfolders)
    
    
    % sub folder path
    sfp = [folder_path '\' subfolders{isf}];
    if contains(sfp, 'gen')
        sfp = [{[sfp '\d01']} {[sfp '\d02']}];
    else
        sfp = {[sfp '\d01']};
    end
    
    
    
    % pdf paths
    for isfp = 1:length(sfp)
        image_path = get_file_paths_targeted(sfp{isfp}, 'intermediate_results.mat')
    
        % subplot
        iplot = iplot+1;
        subplot(1,8,iplot)
        
        % plot
        try
            
            image_path_slashes = strfind(image_path{1}, '\');
            image_path_day = strfind(image_path{1}, '\d0');
            image_path_slashes = max(image_path_slashes(image_path_slashes<image_path_day));
            path_title = image_path{1}(image_path_slashes+1 : image_path_day+3);

            load(image_path{1}, 'initialization')
            imagesc(initialization.neuron.Cn); colormap(bone)
            title(path_title)
            axis([-112.4668  761.4668    0.5000  486.5000]); axis square; axis off
            drawnow
            
        catch
            axis off
            title(path_title)
        end
    end
    
    
    
    
    
    
end
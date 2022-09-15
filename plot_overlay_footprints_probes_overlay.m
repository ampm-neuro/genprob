function plot_overlay_footprints_probes_overlay(folder_path)
% input the folder path containing cnmfe outputs and it will subplot the
% overlay_footprint pdfs from every day, in chronological order

% sub folders
session_chron_reorder_idx = [1 9 15 2 10 16 3 11 17 4 12 18 5 13 19 6 14 20 7 8];
session_numbers = [1:3:19 20];
subfolders = [{'probe01'} {'probe02'} {'probe03'} {'probe04'} {'probe05'} {'probe06'} {'probe07'} {'probe08'}];

% load cell reg file
subj_id = find_subj_id(folder_path);
cell_reg_path = ['C:\Users\ampm1\Documents\MATLAB\generalization\cellreg_data\' subj_id '\cell_reg_' subj_id '.mat'];
load(cell_reg_path, 'cell_registered_struct');

% figure
figure;
          
% iterate through the sub folders          
iplot_idx = 0;
iplot = [1 9 17 25 2 10 18 26 3 11 19 27 4 12 20 28 5 13 21 29 6 14 22 30 7 15 23 31 8 16 24 32];
local_paths = [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6 7 7 7 7 8 8 8 8];
for isf = 1:length(iplot)
    
    % sub folder path
    sfp = [folder_path '\' subfolders{local_paths(isf)}];
    if contains(sfp, 'gen')
        sfp = [{[sfp '\d01']} {[sfp '\d02']}];
    else
        sfp = {[sfp '\d01']};
    end
    
    
    % pdf paths
    for isfp = 1:length(sfp)
    
        % subplot
        iplot_idx = iplot_idx+1;
        subplot(4,8,iplot(iplot_idx)); hold on

        % Images for first two rows
        if iplot(iplot_idx) < 17
        
            % image path
            image_path = get_file_paths_targeted(sfp{isfp}, 'intermediate_results.mat');
            image_path_slashes = strfind(image_path{1}, '\');
            image_path_day = strfind(image_path{1}, '\d0');
            image_path_slashes = max(image_path_slashes(image_path_slashes<image_path_day));
            %path_title = image_path{1}(image_path_slashes+1 : image_path_day+3);
            path_title = image_path{1}(image_path_slashes+1 : image_path_day-1);

            % load image
            load(image_path{1}, 'initialization')
            
            % plot image
            imagesc(initialization.neuron.Cn); colormap(bone)
            
            % aesthetics
            title(path_title)
            axis([-112.4668  761.4668    0.5000  486.5000]); axis square; axis off
            set(gca, 'YDir','reverse')
        end
        
        % Raw footprints for second row
        if iplot(iplot_idx) > 8 && iplot(iplot_idx) < 17
        
            % load
            footprints_path = get_file_paths_targeted(sfp{isfp}, 'cnmfe_out.mat');
            load(footprints_path{1}, 'A', 'Cn', 'rev_traces')

            % preallocate
            footprints = nan(sum(rev_traces), size(Cn,1), size(Cn,2));

            % approved traces
            apr_trc = find(rev_traces==1);

            % reshape into matrices
            for itrace = 1:length(apr_trc)
                % current approved trace
                cat = apr_trc(itrace);
                % reshape
                footprints(itrace,:,:) = reshape(A(:,cat), size(Cn,1), size(Cn,2));
            end

            % plot footprints
            %trace_footprints(footprints, [186 032 013]./255)
            trace_footprints(footprints, [0 0 1])
            
            % aesthetics
            title(path_title)
            axis([-112.4668  761.4668    0.5000  486.5000]); axis square; axis off
            set(gca, 'YDir','reverse')
        end
        
        % Footprints only for third row
        if iplot(iplot_idx) >= 17 && iplot(iplot_idx) < 25

            % adjacent session numbers
            session_num = session_chron_reorder_idx(session_numbers([local_paths(isf)]))
            
            % plot
            trace_footprints_sessions(cell_registered_struct.spatial_footprints_corrected(session_num), [0 0 1])
            legend off
            
            % aesthetics
            title(num2str(session_num))
            axis([-112.4668  761.4668    0.5000  486.5000]); axis square; axis off
            set(gca, 'YDir','reverse')
        end
        
        % Registered, overlapping footprints for fourth row
        if iplot(iplot_idx) >= 25 && iplot(iplot_idx) < 32

            % adjacent session numbers
            session_comparison = session_chron_reorder_idx(session_numbers([local_paths(isf) local_paths(isf)+1]))
            
            % plot
            trace_footprints_sessions_overlap(cell_registered_struct.spatial_footprints_corrected(session_comparison), cell_registered_struct.cell_to_index_map(:,session_comparison), 0)
            legend off
            
            % aesthetics
            title(num2str(session_comparison))
            axis([-112.4668  761.4668    0.5000  486.5000]); axis square; axis off
            set(gca, 'YDir','reverse')
        end
        
        % last box empty
        if iplot(iplot_idx)==32
            axis off
        end
    end
    
drawnow    
end
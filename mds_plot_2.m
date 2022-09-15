function [mds_coords, stress] = mds_plot_2(mtx, dims, plot_idx, downsamp)
% plots output of mds dim reduction
% dims indicates the number of dimensions after reduction, try 2


% downsample
if size(mtx,2) < size(plot_idx,2)
    plot_idx = plot_idx(1:downsamp:end);
else
    mtx = mtx(1:downsamp:end, 1:downsamp:end);
    plot_idx = plot_idx(1:downsamp:end);
end
upi = unique(plot_idx)';

% colors for plotting
if exist('plot_idx', 'var')
    
    % main colors (1 per session)
    colors_temp = hsv(length(unique(plot_idx))+1);
    colors = nan(size(colors_temp,1)-1,3);
    for icolor = 1:size(colors,1)
        colors(icolor,:) = mean(colors_temp(icolor:icolor+1,:));
    end
    
    % fading between main colors (for lines)
    colors_2 = nan(size(mtx,1)-1,3);
    for icolor = 1:length(upi)
    
        if icolor == 1
            color_hold = [linspace(colors_temp(icolor,1), colors_temp(icolor+1,1), sum(plot_idx==upi(icolor)))' ...
                          linspace(colors_temp(icolor,2), colors_temp(icolor+1,2), sum(plot_idx==upi(icolor)))' ...
                          linspace(colors_temp(icolor,3), colors_temp(icolor+1,3), sum(plot_idx==upi(icolor)))'];
            colors_2(plot_idx==upi(icolor),:) = color_hold;
        else
            color_hold = [linspace(colors_temp(icolor,1), colors_temp(icolor+1,1), sum(plot_idx==upi(icolor))+1)' ...
                          linspace(colors_temp(icolor,2), colors_temp(icolor+1,2), sum(plot_idx==upi(icolor))+1)' ...
                          linspace(colors_temp(icolor,3), colors_temp(icolor+1,3), sum(plot_idx==upi(icolor))+1)'];
            colors_2(plot_idx==upi(icolor),:) = color_hold(2:end, :);
        end
    end
    line_ct = 0;
end



% computes mds coords
[mds_coords, stress] = mdscale(mtx, dims);

% plots first 3 coords
%
if dims>=3
    figure; hold on;
    legend_trick(colors, '.')
    
    if exist('plot_idx', 'var')
        pi_means = nan(length(upi),dims);
        for iidx = 1:length(upi)
              
            % current session data
            cpi = plot_idx==upi(iidx);
            mds_coords_local = mds_coords(cpi,:);
            if iidx<length(upi)
                colors_2_local = colors_2(cpi,:);
            else
                colors_2_local = colors_2(cpi(1:end-1),:);
            end
            
            % plot lines
            for idot = 1:size(mds_coords_local,1)-1
                line_ct = line_ct+1;
                line_pts = idot : idot+1;
                %colors_2(line_ct,:)
                plot3(mds_coords_local(line_pts,1), mds_coords_local(line_pts,2), mds_coords_local(line_pts,3), '-', 'color', colors_2(line_ct,:))
            end
            
            % plot dots
            plot3(mds_coords(cpi,1), mds_coords(cpi,2), mds_coords(cpi,3), '.', 'color', colors(iidx,:), 'markersize', 15)

            % current session means (dots and lines)
            pi_means(iidx,:) = mean(mds_coords_local,1);
            if ismember(iidx, [1:3:19 20])
                plot3(pi_means(iidx,1), pi_means(iidx,2), pi_means(iidx,3), '.', 'color', colors(iidx,:), 'markersize', 50)
            else
                plot3(pi_means(iidx,1), pi_means(iidx,2), pi_means(iidx,3), '.', 'color', colors(iidx,:), 'markersize', 50)   
            end
            if iidx>1
                plot3(pi_means(iidx-1:iidx,1), pi_means(iidx-1:iidx,2), pi_means(iidx-1:iidx,3), 'k-', 'linewidth', 3)
            end
        end
    else
        plot3(mds_coords(:,1), mds_coords(:,2), mds_coords(:,3), 'o')
    end
    axis equal
    axis square
    set(gca,'TickLength',[0, 0]); box off;
    
    
    % legend
    legend_string = cell(1, length(upi));
    for ilegend_str = 1:length(legend_string)
       legend_string{ilegend_str} = ['Session ' num2str(upi(ilegend_str))];
    end
    legend(legend_string)
end
%}

% plots first 2 coords

if dims==2
    
    f1 = figure; hold on;
    legend_trick(colors, '.')
    
    if exist('plot_idx', 'var')
        pi_means = nan(length(upi),dims);
        for iidx = 1:length(upi)
              
            % current session data
            cpi = plot_idx==upi(iidx);
            mds_coords_local = mds_coords(cpi,:);
            if iidx<length(upi)
                colors_2_local = colors_2(cpi,:);
            else
                colors_2_local = colors_2(cpi(1:end-1),:);
            end
            
            % plot lines
            for idot = 1:size(mds_coords_local,1)-1
                line_ct = line_ct+1;
                line_pts = idot : idot+1;
                %colors_2(line_ct,:)
                figure(f1)
                plot(mds_coords_local(line_pts,1), mds_coords_local(line_pts,2), '-', 'color', colors_2(line_ct,:))

            end
            
            % plot dots
            plot(mds_coords(cpi,1), mds_coords(cpi,2), '.', 'color', colors(iidx,:), 'markersize', 15)

            % current session means (dots and lines)
            pi_means(iidx,:) = mean(mds_coords_local,1);
            if iidx>1
                plot(pi_means(iidx-1:iidx,1), pi_means(iidx-1:iidx,2), 'k-', 'linewidth', 3)
            end
            plot(pi_means(iidx,1), pi_means(iidx,2), '.', 'color', colors(iidx,:), 'markersize', 50)
        end
    else
        plot(mds_coords(:,1), mds_coords(:,2), 'o')
    end
    axis equal
    axis square
    set(gca,'TickLength',[0, 0]); box off;
    
    
    % legend
    legend_string = cell(1, length(upi));
    for ilegend_str = 1:length(legend_string)
       legend_string{ilegend_str} = ['Session ' num2str(upi(ilegend_str))];
    end
    legend(legend_string)
    legend('location', 'northeastoutside')
    
    
    
end
    
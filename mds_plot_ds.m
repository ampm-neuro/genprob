function [mds_coords, stress] = mds_plot_ds(mtx, dims, downsample, plot_idx)
% plots output of mds dim reduction
% dims indicates the number of dimensions after reduction, try 2

if exist('plot_idx', 'var')
    colors = parula(length(unique(plot_idx)));
end


% downsample
mtx = mtx(1:downsample:end, :);
plot_idx = plot_idx(1:downsample:end);
upi = unique(plot_idx)';


% computes pairwise distances
distances_mtx = distance_window(mtx', mtx', 2);
%distances_mtx(isinf(distances_mtx)) = median(distances_mtx(~isinf(distances_mtx)));
distances_mtx(isinf(distances_mtx)) = nan;
distances_mtx = inpaint_nans(distances_mtx);
distances_mtx(abs(distances_mtx-distances_mtx') > 10*eps*max(distances_mtx(:))) = 0.0001;

% computes mds coords
[mds_coords, stress] = mdscale(distances_mtx, dims);

% plots first 3 coords
%{
if dims>=3
    figure; hold on;
    legend_trick(colors, 'o')
    plot3(mds_coords(:,1), mds_coords(:,2), mds_coords(:,3), 'k-')
    if exist('plot_idx', 'var')
        for iidx = 1:length(upi)
            cpi = plot_idx==plot_idx(iidx);
            plot3(mds_coords(cpi,1), mds_coords(cpi,2), mds_coords(cpi,3), 'o', 'color', colors(iidx,:))
        end
    else
        plot3(mds_coords(:,1), mds_coords(:,2), mds_coords(:,3), 'o')
    end
    axis equal
    axis square
    set(gca,'TickLength',[0, 0]); box off;
end
%}

% plots first 2 coords
if dims==2
    figure; hold on;
    
    if exist('plot_idx', 'var')
        pi_means = nan(length(upi),dims);
        for iidx = 1:length(upi)
            cpi = plot_idx==upi(iidx);
            plot(mds_coords(cpi,1), mds_coords(cpi,2), 'o', 'color', colors(iidx,:))
            pi_means(iidx,:) = mean(mds_coords(cpi,:),1);
            if iidx>1
                plot(pi_means(iidx-1:iidx,1), pi_means(iidx-1:iidx,2), 'k-', 'linewidth', 3)
            end
            plot(mean(mds_coords(cpi,1)), mean(mds_coords(cpi,2)), '.', 'color', colors(iidx,:), 'markersize', 50)
        end
    else
        plot(mds_coords(:,1), mds_coords(:,2), 'k-')
        plot(mds_coords(:,1), mds_coords(:,2), 'o')
    end
    axis equal
    axis square
    set(gca,'TickLength',[0, 0]); box off;
end
    
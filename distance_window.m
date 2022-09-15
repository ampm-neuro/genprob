function distances = distance_window(starting_points, ending_points, plot_line_or_mtx)
% Computes and plots euclidean distance from each COLUMN in 'start_points'  
% to each COLUMN in 'ending_points'. 
% 
% 'distances' is a matrix with 1 row per column in 'ending_points'
%
% plot 0 does not make a plot
% plot 1 shows 1 line for each row in 'distances'
% plot 2 makes a pairwise heatmap

% downsample
%starting_points = starting_points(:, 1:50:end);
%ending_points = ending_points(:, 1:50:end);


% colors
colors = parula(size(ending_points,2));

% preallocate distances
distances = nan(size(ending_points,2), size(starting_points,2));

% iterate through ending points
for iep = 1:size(ending_points,2)
    
    % iterate through starting points
    for isp = 1:size(starting_points,2)
       
        % points
        start_pt = starting_points(:, isp);
        end_pt = ending_points(:, iep);

        % dimensionaltiy correction
        dim_correct = sqrt(sum(start_pt>0 & end_pt>0));
        
        % compute and load distance
        distances(iep, isp) = sqrt(sum((start_pt - end_pt) .^ 2))/dim_correct;
        
    end
end

%% plot


if plot_line_or_mtx==1
    figure; hold on

legend_trick(colors, '-')

for iep = 1:size(ending_points,2)
    plot(distances(iep,:), '-', 'color', colors(iep,:))
end

% aesthetics
xlabel('window')
ylabel('distance')
set(gca,'TickLength',[0, 0]); box off;
xlim([0.5 size(distances,2)+0.5])

% legend
legend_strings = [];
for isesh = 1:size(ending_points,2)
   legend_strings = [legend_strings; {num2str(session_input(isesh))}];
end
legend(legend_strings)


elseif plot_line_or_mtx==2
    
    figure; imagesc(distances)
    
    if size(distances,1)==size(distances,2)
        axis square
    end
    set(gca,'TickLength',[0, 0]); box off;
    colorbar
    
    
end
    



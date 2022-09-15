function [all_dispersion, all_reliability] = all_cells_infoContent(neurons, session_mtx_cell)
% compute and plot average infocontent

all_dispersion = [];
all_reliability = [];

%
for icell = neurons

    all_dispersion = [all_dispersion; ...
        one_cell_infoContent(icell, session_mtx_cell)];
    all_reliability = [all_reliability; ...
        one_cell_reliability(icell, session_mtx_cell)'];
end
%}


%{
% errorbar dispersion
eb_in = cell(1, size(all_dispersion,2));
for isesh = 1:size(all_dispersion,2)
    eb_in{isesh} = all_dispersion(:, isesh);
end
figure; hold on
errorbar_plot(eb_in)
ylabel('Dispersion')


% errorbar reliability
eb_in = cell(1, size(all_reliability,2));
for isesh = 1:size(all_reliability,2)
    eb_in{isesh} = all_reliability(:, isesh);
end
figure; hold on
errorbar_plot(eb_in)
ylabel('Reliability')
%}

%% figure
%{
% can be copy pasted into command line
num_time_bins = size(session_mtx_cell{1}{1},2);
divisor = 10; %set in dispersion calculations (function: contrast_score)
colors = parula(20); 
figure; hold on; 
legend_trick(colors, '.')

sesh_means = nan(20,2); 
for isesh = 1:size(all_reliability,2) 
    sesh_means(isesh,:) = nanmean([all_reliability(:,isesh), 1-all_dispersion(:,isesh)./ceil(num_time_bins/divisor)]); 
    plot(all_reliability(:,isesh), 1-all_dispersion(:,isesh)./(num_time_bins/divisor), 'o', 'color', colors(isesh,:)); 
end

for isesh = 1:size(all_reliability,2)
    plot(sesh_means(isesh,1), sesh_means(isesh,2), '.', 'markersize', 60, 'color', colors(isesh,:));
    
    if isesh>1
        plot(sesh_means(1:isesh,1), sesh_means(1:isesh,2), 'k-', 'linewidth', 5);
    end
end

axis([0 1 0 1])
axis square
set(gca,'TickLength',[0, 0]); box off;
xticks([0:.2:1])
yticks([0:.2:1])

xlabel('Reliability')
ylabel('Specificity')

% fit line
hold on; fit_line(all_reliability(:), 1-all_dispersion(:)./(num_time_bins/divisor), 0)
legend_str = [];
for i = 1:20
    legend_str = [legend_str, {num2str(i)}];
end
legend(legend_str)
legend('location', 'northeastoutside')
%}



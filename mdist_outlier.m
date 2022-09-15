function outlier_index = mdist_outlier(xy_vects, mdistance_threshold)
% computes mdist from each point to the others
% outputs index of points > mdistance_threshold
% input is a two column matrix, column 1 is x positions, c2 is y

% compute distances
mdists = nan(size(xy_vects,1),1);
for ipt = 1:size(xy_vects, 1)
    mdists(ipt) = mahal(xy_vects(ipt,:), xy_vects(setdiff(1:size(xy_vects,1), ipt), :))./sqrt(size(xy_vects,2));
end

% outliers
outlier_index = mdists>mdistance_threshold;

disp(['Outliers found: ' num2str(find(outlier_index==1)')])

% plot
figure; hold on
for ipt=1:size(xy_vects,1)
    plot(xy_vects(ipt,1),xy_vects(ipt,2), 'wo'); 
    if outlier_index(ipt)==0
        text(xy_vects(ipt,1),xy_vects(ipt,2), num2str(ipt));
    else
        text(xy_vects(ipt,1),xy_vects(ipt,2), num2str(ipt), 'Color', [1 0 0]);
    end
end

% axis
kept_mean = nanmean(xy_vects(~outlier_index,:));
min_max_x = [min(xy_vects(:,1)) max(xy_vects(:,1))];
    max_error_x = max(abs([kept_mean(1)-min_max_x(1) kept_mean(1)-min_max_x(2)]))*1.15;
min_max_y = [min(xy_vects(:,2)) max(xy_vects(:,2))];
    max_error_y = max(abs([kept_mean(2)-min_max_y(1) kept_mean(2)-min_max_y(2)]))*1.15;
axis([kept_mean(1)-max_error_x kept_mean(1)+max_error_x kept_mean(2)-max_error_y kept_mean(2)+max_error_y])

% aesthetics
set(gca,'TickLength',[0, 0]); box off;
axis square
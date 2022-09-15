function [trl_nums, freqs, wait_times, rich_idx, TTL_idx] = wait_times_sesh_plot(trl_mtx)
%plots the mean wait times at each frequency


% all wait times and corresponding frequencies
trl_nums = (1:size(trl_mtx,1))';
pt_idx = trl_mtx(:,3)==0;
trl_mtx = trl_mtx(pt_idx,:);
trl_nums = trl_nums(pt_idx);
freqs = floor(trl_mtx(:,2));
rich_idx = rich_trl_idx(trl_mtx);
wait_times = trl_mtx(:,12);

if size(trl_mtx,2)>12
    TTL_idx = trl_mtx(:,13);
else
    TTL_idx = zeros(size(trl_mtx(:,12)));
end


% plot
hold on

% colors
colors = [.4 .4 .4; 0.05 0.29 0.65];

% legend prep
plot(-1:-1+realmin, -1:-1+realmin, '.', 'color', colors(1,:), 'markersize', 24);
plot(-1:-1+realmin, -1:-1+realmin, '.', 'color', colors(2,:), 'markersize', 24);
plot(-1:-1+realmin, -1:-1+realmin, 'o', 'color', colors(1,:));
plot(-1:-1+realmin, -1:-1+realmin, 'o', 'color', colors(2,:));


for itrl = 1:size(trl_mtx,1)
    if rich_idx(itrl)==1 & TTL_idx(itrl)==0
        plot(trl_nums(itrl), wait_times(itrl), '.', 'color', colors(1,:), 'markersize', 24)
    elseif rich_idx(itrl)==1 & TTL_idx(itrl)==1
        plot(trl_nums(itrl), wait_times(itrl), '.', 'color', colors(2,:), 'markersize', 24)
    elseif rich_idx(itrl)==0 & TTL_idx(itrl)==0
        plot(trl_nums(itrl), wait_times(itrl), 'o', 'color', colors(1,:))
    elseif rich_idx(itrl)==0 & TTL_idx(itrl)==1
        plot(trl_nums(itrl), wait_times(itrl), 'o', 'color', colors(2,:))
    end
end



% aesthetics
set(gca,'TickLength',[0, 0]);
ylabel('Wait Durations (s)')
ylim_hold = ylim; ylim([0 60]);
title('Waits')
xlim([-2.01 125.01])
xlabel('Trial number')
legend({'RichTone LaserOFF','RichTone LaserON','PoorTone LaserOFF','PoorTone LaserON'}, 'location', 'northeastoutside')
    




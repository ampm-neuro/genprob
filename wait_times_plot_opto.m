function wait_times_plot_opto(trl_mtx, varargin)
%plots the mean wait times at each frequency

if nargin == 2
    plot_type = varargin{1};
    current_colors = [.4 .4 .4; 0.05 0.29 0.65; 0.9882 0.7294 0.0118];
elseif nargin == 3
    plot_type = varargin{1};
    current_colors = varargin{2};
else
    plot_type = 1;
    current_colors = distinguishable_colors(2);
end

hold on
legend_trick(current_colors, 'o')
wait_times_plot(trl_mtx(trl_mtx(:,13)==0,:), plot_type, current_colors(1,:));
wait_times_plot(trl_mtx(trl_mtx(:,13)==1,:), plot_type, current_colors(2,:));
wait_times_plot(trl_mtx(trl_mtx(:,13)==2,:), plot_type, current_colors(3,:));
ylim([0 60])

legend({'laser OFF', 'laser ON', 'laser ON'})






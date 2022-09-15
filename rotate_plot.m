function out_points = rotate_plot(deg_ccw, in_points)
% rotate counter clockwise
% orig_pts must be a column of x pts and a col of y pts

% center on origin
origin_correction = -nanmean(in_points,2);
origin_points = in_points + origin_correction;

% rotation matrix
rot_mtx = [cosd(deg_ccw) -sind(deg_ccw); sind(deg_ccw) cosd(deg_ccw)];

% rotate about origin
rot_points = rot_mtx*origin_points;

% return to original position relative to origin
out_points = rot_points - origin_correction;

figure; hold on
plot(in_points(1,:), in_points(2,:))
plot(origin_points(1,:), origin_points(2,:))
plot(rot_points(1,:), rot_points(2,:))
plot(out_points(1,:), out_points(2,:))
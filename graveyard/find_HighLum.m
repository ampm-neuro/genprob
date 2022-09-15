function HighLum_idx = find_HighLum(lum_mean, videoObject, min_sec)
% find LED-on instances in a minimic imaging file

% correct brightness
lum_bounds = [50 64];
lum_bounds_idx = lum_mean > lum_bounds(1) & lum_mean < lum_bounds(2);

% correct duration


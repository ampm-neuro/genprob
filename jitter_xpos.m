function xpos_out = jitter_xpos(xpos, ypos, varargin)
%jitters x position of dotplots

%columnize
ypos = ypos(:);
xpos = xpos(:);

%fixed multiplier input
if nargin == 3
    fixed_mult = varargin{1};
else
    fixed_mult = 1;
end

% distribution stats
num_points = length(ypos);
mean_y = nanmean(ypos);
std_y = nanstd(ypos);
dists_from_mean = abs(ypos - mean_y);
std_from_mean = dists_from_mean/std_y; 
if isnan(std_from_mean); std_from_mean=0; end

% jitter
base_jitter_multiplier = num_points*0.0006 + 0.25;
base_jitter = (rand(length(ypos),1)-0.5).*base_jitter_multiplier;
bulb_correction = std_from_mean./7.5 + (rand(length(ypos),1)-0.5)*0.4;
jitter_update = base_jitter.*(1-bulb_correction);

% output 
xpos_out = xpos + jitter_update.*fixed_mult;
                

function trace_footprints(footprint_mtx, varargin)
% plots outline of all footprints
% footprints are a 3d matrix of footprints from a single session
% cell, ~length, ~width

% input
if ~isempty(varargin)
    plot_color = varargin{1};
else
    plot_color = [];
end

hold on

% iterate through each cell
for icell = 1:size(footprint_mtx,1)
    trace_footprint(squeeze(footprint_mtx(icell,:,:)), plot_color);
end



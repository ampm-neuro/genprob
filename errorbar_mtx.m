function errorbar_mtx(mtx_in, varargin)
% inputs a matrix into errorbar function by translating columns into cells

if size(mtx_in,3)>1
    eb_colors = distinguishable_colors(size(mtx_in,3));
    colors = eb_colors./2;
else
    eb_colors = rand(1,3);
    colors = eb_colors./2;
end

if ~isempty(varargin)
    eb_colors = varargin{1};
    colors = varargin{2};
end

% convert mtx to cell
ebp_cell = cell(size(mtx_in,3), size(mtx_in,2));
for ipage = 1:size(mtx_in,3)
    for icol = 1:size(mtx_in,2)
        ebp_cell{ipage, icol} = mtx_in(:, icol, ipage);
    end
end

hold on
for igrp = 1:size(ebp_cell,1)
    errorbar_plot(ebp_cell(igrp,:), 1, [], colors(igrp,:), eb_colors(igrp,:).*0.75)
end

% one overall errorbar
%errorbar_plot_multi(ebp_cell)

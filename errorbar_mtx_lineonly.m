function errorbar_mtx_lineonly(mtx_in, varargin)
% inputs a matrix into errorbar function by translating columns into cells

if size(mtx_in,3)>1
    colors = distinguishable_colors(size(mtx_in,3));
else
    colors = .7.*[1 1 1];
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
    errorbar_plot_lineonly(ebp_cell(igrp,:), [], [], colors(igrp,:), colors(igrp,:).*0.75)
end

% one overall errorbar
%errorbar_plot_multi(ebp_cell)

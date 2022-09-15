function legend_trick(colors, lineshape_str)
% adds small markings to plot that bias legend items

hold on
for icolor = 1:size(colors,1)
    
   plot(-1:-1+realmin, -1:-1+realmin, lineshape_str, 'color', colors(icolor,:)); 
    
end
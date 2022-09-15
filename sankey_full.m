function sankey_full(cell_reg_mtx, sessions)
% makes a sankey diagram showing which cells are active during which
% sessions
% sessions corresponds to columns in cell_reg_mtx

% cell reg active or inactive
cell_reg_mtx = double(cell_reg_mtx(:,sessions)>0);

% number of probes
num_sessions = length(sessions);

%% data of each block
unq_hist = cell(1, num_sessions);
number_cells_in_each_block = cell(1,num_sessions);
number_of_blocks = 0;
for isesh = 1:num_sessions
    
    % relevant history sorted by current session activity
    rel_hist = cell_reg_mtx(:,1:isesh);
    rel_hist(rel_hist==0)=nan;
    %rel_hist = sortrows(rel_hist, [isesh 1:isesh-1]); 
    rel_hist = sortrows(rel_hist, isesh:-1:1); 
    rel_hist(isnan(rel_hist))=0;
    
    % unq histories hard code
    %  
        % preallocate
        unq_hist{isesh} = dec2bin(2^isesh-1:-1:0)-'0';
        unq_hist{isesh} = flipud(sortrows(unq_hist{isesh}, isesh:-1:1)); 
    
    % unique histories
    %unique(rel_hist, 'rows', 'stable')
    %unq_hist{isesh} = unique(rel_hist, 'rows', 'stable');
    num_unq_hist = size(unq_hist{isesh},1);

    % size of each history block
    number_cells_in_each_block{isesh} = nan(1, num_unq_hist); 
    for ihist = num_unq_hist:-1:1
        number_cells_in_each_block{isesh}(ihist) = sum(ismember(rel_hist, unq_hist{isesh}(ihist,:), 'rows'));
    end

    % keep track of number of blocks
    number_of_blocks = number_of_blocks + length(number_cells_in_each_block{isesh});
end

% data of the left panel in each block
data1 = number_cells_in_each_block(1:end-1);

% data of the right panel in each block
data2 = number_cells_in_each_block(2:end);



%% data of flows in each block
data = cell(1, num_sessions-1);
for isesh = 1:num_sessions-1
    
    % preallocate
    data{isesh} = nan(length(number_cells_in_each_block{isesh}), length(number_cells_in_each_block{isesh+1}));
    
    %iterate through unique left panel histories
    rel_hist_lp = cell_reg_mtx(:,1:isesh);
    for ilp = 1:size(data{isesh},1)
        
        unq_hist_lp = unq_hist{isesh};
        idx_lp = ismember(rel_hist_lp, unq_hist_lp(ilp,:), 'rows');
       
        % iterate through unique right panel histories
        rel_hist_rp = cell_reg_mtx(:,1:isesh+1);
        for irp = 1:size(data{isesh},2)
            
            unq_hist_rp = unq_hist{isesh+1};
            idx_rp = ismember(rel_hist_rp, unq_hist_rp(irp,:), 'rows');
            
            % load counts
            data{isesh}(ilp,irp) = sum(idx_lp & idx_rp);
            
        end
    end
end

%% colors

% one color for
% still inactive pool
% each context a cell first becomes active in

%{
all_colors = distinguishable_colors(number_of_blocks);
barcolors = cell(1,length(number_cells_in_each_block));
color_ct = 1;
for isesh = 1:length(number_cells_in_each_block)   
    barcolors{isesh} = all_colors(color_ct:color_ct + length(number_cells_in_each_block{isesh}) -1, :);
    color_ct = color_ct + length(number_cells_in_each_block{isesh});
end
%}

barcolors = cell(1,length(number_cells_in_each_block));
for isesh = 1:length(number_cells_in_each_block) 
    
    num_patches = length(number_cells_in_each_block{isesh});
    half_patches = num_patches/2;
    
    summer_hold = summer(half_patches*2);
    summer_hold = summer_hold(1:half_patches, :);
    summer_hold = flipud(summer_hold);
    winter_hold = bone(half_patches*3);
    winter_hold = winter_hold(half_patches+1:half_patches*2, :);
    winter_hold = flipud(winter_hold);
        
    barcolors{isesh} = [summer_hold; winter_hold];
end



%% plotting


% x-axis
X= 1 : length(number_cells_in_each_block{isesh});


% flow color
c = [.7 .7 .7];

% Panel width
w = 25; 

for j=1:num_sessions-1
    ymax=sankey_yheight(data1{j},data2{j});
    
    if j>1
        %ymax=max(ymax,sankey_yheight(data1{j-1},data2{j-1}));
        y1_category_points=sankey_alluvialflow(data1{j}, data2{j}, data{j}, X(j), X(j+1), y1_category_points,ymax,barcolors{j},barcolors{j+1},w,c);
    else
        y1_category_points=[];
        %ymax=sankey_yheight(data1{j},data2{j});
        y1_category_points=sankey_alluvialflow(data1{j}, data2{j}, data{j}, X(j), X(j+1), y1_category_points,ymax,barcolors{j},barcolors{j+1},w,c);
    end
end

% axis
%{
axes('Color','none','YColor','none');
set(gca,'TickLength',[0, 0]); box off;
xticks(linspace(0,1,length(number_cells_in_each_block)))
xticklabels(sessions)
h=get(gca);
set(h,'xcolor','none')
h.XAxis.Label.Color=[0 0 0];
h.XAxis.Label.Visible='on';
%}
        
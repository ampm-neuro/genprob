function [number_newActive_cells, number_oldActive_cells, number_Active_cells] = sankey_full_probeColor(cell_reg_mtx, sessions)
% makes a sankey diagram showing which cells are active during which
% sessions
% sessions corresponds to columns in cell_reg_mtx

% cell reg active or inactive
cell_reg_mtx = double(cell_reg_mtx(:,sessions)>0);

% number of probes
num_sessions = length(sessions);

%% data of each block
unq_hist = cell(1, num_sessions);
number_newActive_cells = nan(1, num_sessions);
number_oldActive_cells = nan(1, num_sessions);
number_Active_cells = nan(1, num_sessions);
number_cells_in_each_block = cell(1,num_sessions);
number_of_blocks = 0;
for isesh = 1:num_sessions
    
    % relevant history sorted by current session activity
    rel_hist = cell_reg_mtx(:,1:isesh);
    rel_hist(rel_hist==0)=nan;
    %rel_hist = sortrows(rel_hist, [isesh 1:isesh-1]); 
    %rel_hist = sortrows(rel_hist, isesh:-1:1); 
    rel_hist = sortrows(rel_hist, 1:1:isesh);
    rel_hist(isnan(rel_hist))=0;
    
    % unq histories hard code
    %  
        % preallocate
        unq_hist{isesh} = dec2bin(2^isesh-1:-1:0)-'0';
        %unq_hist{isesh} = flipud(sortrows(unq_hist{isesh}, isesh:-1:1));
        %flipud(sortrows(unq_hist{isesh}, isesh:-1:1))
        %flipud(sortrows(unq_hist{isesh}, 1:1:isesh))
        unq_hist{isesh} = flipud(sortrows(unq_hist{isesh}, [isesh 1:isesh-1]));
    
    % unique histories
    %unique(rel_hist, 'rows', 'stable')
    %unq_hist{isesh} = unique(rel_hist, 'rows', 'stable');
    num_unq_hist = size(unq_hist{isesh},1);

    % size of each history block
    number_cells_in_each_block{isesh} = nan(1, num_unq_hist); 
    for ihist = num_unq_hist:-1:1
        number_cells_in_each_block{isesh}(ihist) = sum(ismember(rel_hist, unq_hist{isesh}(ihist,:), 'rows'));
    end
    
    number_newActive_cells(isesh) = sum(ismember(rel_hist, [zeros(1,isesh-1) 1], 'rows'));
    if isesh>1
        number_oldActive_cells(isesh) = sum(rel_hist(:,end)==1 & sum(rel_hist(:,1:end-1),2)>0);
    end
    number_Active_cells(isesh) = sum(rel_hist(:,end)==1);

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

% one color per block, each block has history defined by
% unq_hist{isesh}

%{
all_bcolors = [097 191 050; ...
               186 041 041; ...
               051 153 255; ...
               255 128 000; ...
               143 020 181; ...
               249 221 037; ...
               204 000 102; ...
               120 120 120];
%}
           
all_bcolors = [.3.*[255 255 255];...
    186 41 41; ...
    232 116 0; ...
    232 204 37; ...
    28 142 22; ...
    22 178 247; ...
    48 83 193; ....
    143 020 181];
           
           
all_bcolors = all_bcolors./255;
               
               
               
               
               
    
%distinguishable_colors(8)];
inactive_color = [0 51 102]./255;

barcolors = cell(1,length(sessions));
for isesh = 1:length(sessions) 
    barcolors{isesh} = nan(length(number_cells_in_each_block{isesh}), 3);
    
    for ipatch = 1:length(number_cells_in_each_block{isesh})

        if any(unq_hist{isesh}(ipatch,:)==1)
            barcolors{isesh}(ipatch,:) = all_bcolors(find(unq_hist{isesh}(ipatch,:)==1, 1, 'first'),:);
        else
            barcolors{isesh}(ipatch,:) = inactive_color;
        end
    end
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
        
%% bar charts

figure; bar(number_newActive_cells); title('number new active'); set(gca,'TickLength',[0, 0]); box off;
figure; bar(number_oldActive_cells); title('number old active'); set(gca,'TickLength',[0, 0]); box off;
prop_new = number_newActive_cells./number_Active_cells;
prop_old = number_oldActive_cells./number_Active_cells;
figure; bar([prop_new; prop_old]'); ylim([0 1]); title('proportion new and old active'); set(gca,'TickLength',[0, 0]); box off;

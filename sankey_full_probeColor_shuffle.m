function [number_newActive_cells_shuf, number_oldActive_cells_shuf, number_Active_cells_shuf] = sankey_full_probeColor_shuffle(cell_reg_mtx, sessions, num_shufs)
% makes a sankey diagram showing which cells are active during which
% sessions
% sessions corresponds to columns in cell_reg_mtx

% cell reg active or inactive
cell_reg_mtx = double(cell_reg_mtx(:,sessions)>0);

% number of probes
num_sessions = length(sessions);

%% data of each block
unq_hist = cell(1, num_sessions);
number_newActive_cells_shuf = nan(num_shufs, num_sessions);
number_oldActive_cells_shuf = nan(num_shufs, num_sessions);
number_Active_cells_shuf = nan(num_shufs, num_sessions);
number_cells_in_each_block = cell(1,num_sessions);
number_of_blocks = 0;

% iterate through shuffles
for ishuf = 1:num_shufs
    
    %shuffle
    cell_reg_mtx_shuf = cell_reg_mtx(:,randperm(size(cell_reg_mtx,2)));
    
    %iterate through sessions
    for isesh = 1:num_sessions

        % relevant history sorted by current session activity
        rel_hist = cell_reg_mtx_shuf(:,1:isesh);
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
        if isesh>1
        number_newActive_cells_shuf(ishuf, isesh) = sum(ismember(rel_hist, [zeros(1,isesh-1) 1], 'rows'));    
        number_oldActive_cells_shuf(ishuf, isesh) = sum(rel_hist(:,end)==1 & sum(rel_hist(:,1:end-1),2)>0);
        end
        number_Active_cells_shuf(ishuf, isesh) = sum(rel_hist(:,end)==1);

        % keep track of number of blocks
        number_of_blocks = number_of_blocks + length(number_cells_in_each_block{isesh});
    end
end


%% shuffle plots

% true values
figure; 
[number_newActive_cells, number_oldActive_cells, number_Active_cells] = sankey_full_probeColor(cell_reg_mtx, sessions);

% histograms showing difference between early and late probes
early_probes = [2 3 4];
late_probes = [5 6 7];

% number of new cells
true_val = (sum(number_newActive_cells(early_probes)) - sum(number_newActive_cells(late_probes)));
shuf_early = sum(number_newActive_cells_shuf(:,early_probes),2);
shuf_late = sum(number_newActive_cells_shuf(:,late_probes),2);
shuf_difference = sort(shuf_early-shuf_late);
pval = sum(shuf_difference>true_val)/length(shuf_difference);
figure; hold on; 
histogram(shuf_difference)
plot(true_val.*[1 1], ylim, 'r-')
plot(0.*[1 1], ylim, 'k--')
set(gca,'TickLength',[0, 0]); box off;
ylabel('Number of cells')
xlabel('Number of new cells, early probes - late probes')
title(['number of new cells; pval = ' num2str(pval)])

% number of old cells
true_val = (sum(number_oldActive_cells(early_probes)) - sum(number_oldActive_cells(late_probes)));
shuf_early = sum(number_oldActive_cells_shuf(:,early_probes),2);
shuf_late = sum(number_oldActive_cells_shuf(:,late_probes),2);
shuf_difference = sort(shuf_early-shuf_late);
pval = sum(shuf_difference>true_val)/length(shuf_difference);
figure; hold on; 
histogram(shuf_difference)
plot(true_val.*[1 1], ylim, 'r-')
plot(0.*[1 1], ylim, 'k--')
set(gca,'TickLength',[0, 0]); box off;
ylabel('Number of cells')
xlabel('Number of old cells, early probes - late probes')
title(['number of old cells; pval = ' num2str(pval)])

% proportion of new cells
prop_new_cells = number_newActive_cells./number_Active_cells;
prop_new_cells_shuf = number_newActive_cells_shuf./number_Active_cells_shuf;
true_val = (mean(prop_new_cells(early_probes)) - mean(prop_new_cells(late_probes)));
shuf_early = mean(prop_new_cells_shuf(:,early_probes),2);
shuf_late = mean(prop_new_cells_shuf(:,late_probes),2);
shuf_difference = sort(shuf_early-shuf_late);
pval = sum(shuf_difference>true_val)/length(shuf_difference);
figure; hold on; 
histogram(shuf_difference, 'normalization', 'probability')
plot(true_val.*[1 1], ylim, 'r-')
plot(0.*[1 1], ylim, 'k--')
xlim([-1 1])
set(gca,'TickLength',[0, 0]); box off;
ylabel('Proportion of cells')
xlabel('Mean proportion new cells, early probes - late probes')
title(['Extremeness of new-cell proportions; pval = ' num2str(pval)])

% proportion of old cells
prop_old_cells = number_oldActive_cells./number_Active_cells;
prop_old_cells_shuf = number_oldActive_cells_shuf./number_Active_cells_shuf;
true_val = (mean(prop_old_cells(early_probes)) - mean(prop_old_cells(late_probes)));
shuf_early = mean(prop_old_cells_shuf(:,early_probes),2);
shuf_late = mean(prop_old_cells_shuf(:,late_probes),2);
shuf_difference = sort(shuf_early-shuf_late);
pval = sum(shuf_difference>true_val)/length(shuf_difference);
figure; hold on; 
histogram(shuf_difference, 'normalization', 'probability')
plot(true_val.*[1 1], ylim, 'r-')
plot(0.*[1 1], ylim, 'k--')
xlim([-1 1])
set(gca,'TickLength',[0, 0]); box off;
ylabel('Proportion of cells')
xlabel('Mean proportion new cells, early probes - late probes')
title(['Extremeness of old-cell proportions; pval = ' num2str(pval)])



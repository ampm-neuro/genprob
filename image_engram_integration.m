function image_engram_integration(crm)
% makes plots related to the appearance and inclusion of cells after each
% session (1 session per column in cell_reg_mtx)


%% colors
num_colors = 20;
colors = hsv(num_colors);

%% sort cell_reg_mtx

% set zeros to nans
crm(crm==0) = nan;

% label cells by origin session
for icell = 1:size(crm,1)
    crm(icell, ~isnan(crm(icell,:))) = find(~isnan(crm(icell,:)),1,'first')-1;
end
% sort by first session
crm_new = nan(size(crm));
ct = 0;
for isesh = 1:size(crm,2)
    crm_new(ct+1:ct+sum(crm(:,isesh)==isesh-1), :) = crm(crm(:,isesh)==isesh-1,:);
    ct = ct+sum(crm(:,isesh)==isesh-1);
end
crm = crm_new;

% sort within a session by reactivation probability
session_pct_react = nan(20,1);
for isesh = 1:size(crm,2)
    
    % just cells that originated during this session
    crm_sesh = crm(crm(:,isesh)==isesh-1,:);
    
    % reorder by number of reactivations
    [~, sort_idx] = sort(sum(crm_sesh==isesh-1,2), 'descend');
    crm(crm(:,isesh)==isesh-1,:) = crm_sesh(sort_idx,:);
    
    % load proportion of cells that are ever reactivated
    session_pct_react(isesh) = sum(sum(crm_sesh==isesh-1,2)>1)/size(crm_sesh,1);
end

figure; hold on
%session_pct_react = smooth(session_pct_react,3);
for isesh = 1:size(crm,2)
    ba = bar(isesh, session_pct_react(isesh), 'FaceColor', 'flat');
    ba.CData = colors(isesh,:);
end
title('proportion of origin cells that are reactivated at least once')
set(gca,'TickLength',[0, 0]); box off;
%error

%% plot cell reg mtx


figure; ampm_pcolor(crm')
yticks([])
ylabel('Sessions')
xlabel('Neurons')
set(gcf, 'Position', [903   -32   455   658])
colormap(colors)
box off


%% plot proportion of neurons that are reactivated
sum_react = nan(size(crm,2));
for isesh1 = 1:size(crm,2)
    for isesh2 = 1:size(crm,2)
        sum_react(isesh2, isesh1) = sum(~isnan(crm(:,isesh1)) & crm(:,isesh1)==isesh2-1 & crm(:,isesh1)~=isesh1-1);
    end
end

figure; 
ba = bar((sum_react./sum(~isnan(crm)))', 'stacked', 'FaceColor', 'flat');
for iba_color = 1:size(crm,2)
    ba(iba_color).CData = colors(iba_color,:);
end

set(gca,'TickLength',[0, 0]); box off;
xlim([0.5 20+0.5])
%set(gcf, 'Position', [903   596   455   211])
ylabel('Reactivated portion')
xlabel('Sessions')



%% plot reactivation probability for each session across all remaining sessions
prob_react = nan(1,size(crm,2));
starter_cells = nan(1,size(crm,2)); 
for isesh = 1:size(crm,2)
    crm_sesh = crm(crm(:,isesh)==isesh-1,:);
    starter_cells(isesh) = size(crm_sesh,1);
    cell_probs = nan(1,size(crm_sesh,1));
    for icell = 1:size(crm_sesh,1)
        cell_probs(icell) = sum(crm_sesh(icell,isesh+1:end)==isesh-1)/(size(crm_sesh,2)-isesh);
    end
    prob_react(isesh) = mean(cell_probs);
end

figure; hold on
for isesh = 1:size(crm,2)
    %middle_yval = median(find(crm(:,isesh)==isesh-1));
    % width = num_new_cells-.1;
    num_new_cells(isesh) = sum(crm(:,isesh)==isesh-1);
    h = bar(isesh, prob_react(isesh));
    set(h, 'FaceColor', colors(isesh,:))
end
set(gca,'TickLength',[0, 0]); box off;
%set(gca, 'XDir','reverse')
%set(gca,'view',[90 -90])
xlim([0.5 20+0.5])
ylim([0 1])
%set(gcf, 'Position', [1332 -32 222 658])
xlabel('Sessions')
ylabel('Reactivation probability')


%{

reactivation_probs = sum_react./repmat(starter_cells', 1, length(crm));
reactivation_probs(reactivation_probs==0) = nan;
figure;
first_last = nan(size(crm,2), 2);
subplot(1,2,1); hold on
for isesh = 1:size(crm,2)
    plot(reactivation_probs(isesh,:), 'color', colors(isesh,:))
    
    if sum(~isnan(reactivation_probs(isesh,:)))<2
        continue
    end
    
    first_sesh = reactivation_probs(isesh,find(~isnan(reactivation_probs(isesh,:)), 1, 'first'));
    last_sesh = reactivation_probs(isesh,find(~isnan(reactivation_probs(isesh,:)), 1, 'last'));
    if ~isempty(first_sesh) && ~isempty(last_sesh)
        first_last(isesh,:) = [first_sesh last_sesh];
    end
end
xlabel('Sessions')
ylabel('Reactivation probability')
set(gca,'TickLength',[0, 0]); box off;
ylim([0 1])

% change over time
subplot(1,2,2)
errorbar_plot([{first_last(:,1)} {first_last(:,2)}], 1, [1 2], colors);
xticks([1 2])
ylim([0 1])
%}

%% plot total number of active cells and number of cells originating in each session
total_active_cells = nan(1, size(crm,2));
for isesh = 1:size(crm,2)
    total_active_cells(isesh) = sum(~isnan(crm(:, isesh)));
end
figure; bar([total_active_cells; starter_cells]')
set(gca,'TickLength',[0, 0]); box off;
xlabel('Session')
ylabel('Number of neurons')
xlim([0.5 20+0.5])



figure; bar([starter_cells; total_active_cells-starter_cells]', 'stacked')

proportions = [starter_cells; total_active_cells-starter_cells]./sum([starter_cells; total_active_cells-starter_cells]);
figure; bar(proportions', 'stacked')

% try as a stacked bar
%
figure; hold on
for isesh1 = 1:size(crm,2)      
    colors_temp = [colors(isesh1,:); flipud(colors(setdiff(1:size(crm,2), isesh1),:))]; 
    single_stack = [starter_cells(isesh1) flipud(sum_react(setdiff(1:size(sum_react,1), isesh1), isesh1))'];
    ba = bar([isesh1 nan], [single_stack; nan(size(single_stack))], 'stacked', 'FaceColor', 'flat');
    for isesh2 = 1:size(crm,2)
        ba(isesh2).CData = colors_temp(isesh2,:);
    end
end


sum_act = nan(size(crm,2));
for isesh1 = 1:size(crm,2)
    for isesh2 = 1:size(crm,2)
        sum_act(isesh2, isesh1) = sum(~isnan(crm(:,isesh1)) & crm(:,isesh1)==isesh2-1);
    end
end
figure; 
ba = bar(sum_act', 1.0,  'stacked', 'FaceColor', 'flat');
for iba_color = 1:size(crm,2)
    ba(iba_color).CData = colors(iba_color,:);
end
set(gca,'TickLength',[0, 0]); box off;
xlabel('Session')
ylabel('Number of neurons')
xlim([0.5 20+0.5])
xticks(1:size(crm,2))
%}



%% plot number of new cells vs their probability of reactivation

figure; hold on
for isesh=1:19
    current_x(isesh) = num_new_cells(isesh)./total_active_cells(isesh);
    current_y(isesh) = session_pct_react(isesh);
    if isesh>1
        plot(current_x(isesh), current_y(isesh), '.', 'color', colors(isesh,:), 'markersize', 50);
        if isesh>2
        plot(current_x(isesh-1 : isesh), current_y(isesh-1 : isesh), '-', 'color', colors(isesh,:))
        end
    end
end 
xlabel('proportion of active cells that are new')
ylabel('proportion new cells that will be reactivated')
axis([0 1 0 1]); axis square
set(gca,'TickLength',[0, 0]); box off;




%% Matrix of proportion of cells common to both sessions
figure;
rm = nan(size(crm,2)); 
for i1 = 1:size(crm,2)
    for i2 = 1:size(crm,2) 
        common_cell_sum = sum(~isnan(crm(:,i1)) & ~isnan(crm(:,i2)));
        uncommon_cell_sum = sum([sum(~isnan(crm(:,i1)) & isnan(crm(:,i2))) sum(isnan(crm(:,i1)) & ~isnan(crm(:,i2)))]);
        rm(i1,i2) = common_cell_sum/(common_cell_sum + uncommon_cell_sum);
    end
end
ampm_pcolor(rm); title('Proportion of cells common to both sessions'); axis square; caxis([0 .5]); colorbar
%figure; imagesc(rm)

% plot off-diaganol to show continuity correlation 
%
rm_minor = rm(2:end, 1:end-1);
figure; plot(rm_minor(logical(eye(size(rm_minor,1)))))
set(gca,'TickLength',[0, 0]); box off;
xlim([0.5 size(rm_minor,1)+0.5])
ylim([0 1])
xlabel('session comparison')
ylabel('proportion of shared neurons')
%}




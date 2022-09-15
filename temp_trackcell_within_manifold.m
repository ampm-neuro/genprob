
%figure; imagesc(norm_mtx(nanmean(trial_activity_mtx_14,3)')')

mtx = nanmean(trial_activity_mtx_11,3); 
mtx(:,isnan(mtx(1,:)))=[]; 
mtx = mtx(:, 1:430);
mtx = norm_mtx(mtx')';
[mtx, sort_idx, peakpos] = sort_rows_by_peak(mtx); 
figure; imagesc(mtx); title(11);
hold on; 
plot(peakpos(cell_regist_mtx(36,11)).*[1 1], ylim, 'g--')
plot(xlim, find(sort_idx==cell_regist_mtx(36,11)).*[1 1], 'g--')
cumsumtses = cumsum(tses).*100;
for icst = 1:3
    plot(cumsumtses(icst).*[1 1], ylim, 'r-')
end
set(gca,'TickLength',[0, 0]); box off;axis square


mtx = nanmean(trial_activity_mtx_12,3); 
mtx(:,isnan(mtx(1,:)))=[]; 
mtx = mtx(:, 1:430);
mtx = norm_mtx(mtx')';
[mtx, sort_idx, peakpos] = sort_rows_by_peak(mtx); 
figure; imagesc(mtx); title(12);
hold on; 
plot(peakpos(cell_regist_mtx(36,12)).*[1 1], ylim, 'g--')
plot(xlim, find(sort_idx==cell_regist_mtx(36,12)).*[1 1], 'g--')
cumsumtses = cumsum(tses).*100;
for icst = 1:3
    plot(cumsumtses(icst).*[1 1], ylim, 'r-')
end
set(gca,'TickLength',[0, 0]); box off;axis square



mtx = nanmean(trial_activity_mtx_13,3); 
mtx(:,isnan(mtx(1,:)))=[]; 
mtx = mtx(:, 1:430);
mtx = norm_mtx(mtx')';
[mtx, sort_idx, peakpos] = sort_rows_by_peak(mtx); 
figure; imagesc(mtx); title(13);
hold on; 
plot(peakpos(cell_regist_mtx(36,13)).*[1 1], ylim, 'g--')
plot(xlim, find(sort_idx==cell_regist_mtx(36,13)).*[1 1], 'g--')
cumsumtses = cumsum(tses).*100;
for icst = 1:3
    plot(cumsumtses(icst).*[1 1], ylim, 'r-')
end
set(gca,'TickLength',[0, 0]); box off;axis square


mtx = nanmean(trial_activity_mtx_14,3); 
mtx(:,isnan(mtx(1,:)))=[]; 
mtx = mtx(:, 1:430);
mtx = norm_mtx(mtx')';
[mtx, sort_idx, peakpos] = sort_rows_by_peak(mtx); 
figure; imagesc(mtx); title(14);
hold on; 
plot(peakpos(cell_regist_mtx(36,14)).*[1 1], ylim, 'g--')
plot(xlim, find(sort_idx==cell_regist_mtx(36,14)).*[1 1], 'g--')
cumsumtses = cumsum(tses).*100;
for icst = 1:3
    plot(cumsumtses(icst).*[1 1], ylim, 'r-')
end
set(gca,'TickLength',[0, 0]); box off;axis square


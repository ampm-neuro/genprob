function com = cell_overlap_mtx(cell_regist_mtx_cell, session_chron_reorder_idx, varargin)
% compute proportion of cells that are common to both sessions
% works across mice, 1 crm mtx per subj


% hivar
% load('cell_reg_all.mat')
% cell_overlap_mtx([{cell_regist_mtx_651049m1} {cell_regist_mtx_658648m2}], session_chron_reorder_idx)
% cell_overlap_mtx([{cell_regist_mtx_683472m2} {cell_regist_mtx_683472m3} ], session_chron_reorder_idx)
% session_chron_reorder_idx = [1 9 15 2 10 16 3 11 17 4 12 18 5 13 19 6 14 20 7 8];

% preallocate sums
common_cell_sum_mtx = zeros(size(cell_regist_mtx_cell{1},2));
uncommon_cell_sum_mtx = zeros(size(cell_regist_mtx_cell{1},2));

% iterate through subjects
for isubj = 1:length(cell_regist_mtx_cell)

    % reorder crm
    try
        [~,scri] = sort(session_chron_reorder_idx);
        [~,scri] = sort(scri);
        cell_regist_mtx_chron = cell_regist_mtx_cell{isubj}(:,scri);
    catch
        cell_regist_mtx_chron = cell_regist_mtx_cell{isubj}(:,setdiff(session_chron_reorder_idx, 20, 'stable'));
    end
    
    % iterate through sessions
    for i1 = 1:size(cell_regist_mtx_chron,2)
        for i2 = 1:size(cell_regist_mtx_chron,2)
            common_cell_sum_mtx(i1,i2) = common_cell_sum_mtx(i1,i2) + sum(cell_regist_mtx_chron(:,i1)>0 & cell_regist_mtx_chron(:,i2)>0);
            uncommon_cell_sum_mtx(i1,i2) = uncommon_cell_sum_mtx(i1,i2) + sum([sum(cell_regist_mtx_chron(:,i1)>0 & cell_regist_mtx_chron(:,i2)==0) sum(cell_regist_mtx_chron(:,i1)==0 & cell_regist_mtx_chron(:,i2)>0)]);
        end
    end

end

% overlap 
com = common_cell_sum_mtx./(common_cell_sum_mtx + uncommon_cell_sum_mtx);

figure;
ampm_pcolor(com); 
if ~isempty(varargin)
    title(['Percentage of cells common to both sessions (' varargin{1} ')'])
else
    title('Percentage of cells common to both sessions'); 
end
axis square; caxis([0 .5]); colorbar
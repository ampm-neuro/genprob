function reliability_corm = ALL_reliability_corm(session_mtx_cell, cell_regist_mtx)
% computes the average trial by trial correlation matrix from all cells
% constrain sessions by limiting which session_mtx_cell cells are input

first_last = 1;

reliability_corm = [];

% iterate through sessions
for isesh = 1:length(session_mtx_cell)
        
    % determine which cells to include
    if first_last==1 % first appearance only
    
        % unique cells
        if isesh==1
            % all active cells
            uc_idx = cell_regist_mtx(:,isesh)>0;
            unq_cells = cell_regist_mtx(uc_idx, isesh);
        else
            % cells that haven't been active before this session
            uc_idx = sum(cell_regist_mtx(:,isesh-1:-1:1),2)==0 & cell_regist_mtx(:,isesh)>0;
            unq_cells = cell_regist_mtx(uc_idx, isesh);
        end
    
    elseif first_last ==2 % last appearance only
    
    
        % unique cells
        if isesh==size(cell_regist_mtx,2)
            uc_idx = cell_regist_mtx(:,isesh)>0;
            unq_cells = cell_regist_mtx(uc_idx, isesh);
        else
            uc_idx = sum(cell_regist_mtx(:, isesh+1:1:size(cell_regist_mtx,2)),2)==0 & cell_regist_mtx(:,isesh)>0;
            unq_cells = cell_regist_mtx(uc_idx, isesh);
        end
    

    else % include redundancies
        error('incomplete code')
    end
    
    % compute correlations
    for icell = unq_cells' % local (non-universal) session cell ids

        % limit trials
        activity = session_mtx_cell{isesh}{icell}(1:41,:);

        % correlations
        act_rel = nan(size(activity,1));
        for itrl1 = 1:size(activity,1)
            for itrl2 = 1:size(activity,1)
                if itrl1>=itrl2
                    continue
                end
                act_rel(itrl1,itrl2) = corr(activity(itrl1,:)', activity(itrl2,:)');
            end
        end

        % load each cell
        reliability_corm = cat(3, reliability_corm, act_rel);
    end
end
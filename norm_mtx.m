function mtx_norm = norm_mtx(mtx)
%norm 0 to 1
mtx = mtx - repmat(min(mtx), size(mtx,1), 1);
mtx_norm = mtx./max(mtx);

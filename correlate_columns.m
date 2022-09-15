function corr_mtx = correlate_columns(mtx)
% pairwise correlation of all columns

% preallcate correlation matrix
corr_mtx = nan(size(mtx,2), size(mtx,2));

% iterate through columns
for icol1 = 1:size(mtx,2)
    for icol2 = 1:size(mtx,2)
        corr_mtx(icol1,icol2) = corr(mtx(:,icol1), mtx(:,icol2));
    end
end
function var_out = variance(mtx)
% var without the unbiased estimator

if size(mtx,1)==1
    mtx = mtx(:);
end

% col means
col_means = mean(mtx);

% subtract by col means
mtx_devs = mtx - repmat(col_means, size(mtx,1),1);

% all sqrs
mtx_devs = mtx_devs.^2;

% col sums
col_dev_sum = sum(mtx_devs,1);

% col dev means
var_out = col_dev_sum./size(mtx,1);

% sqrt root
%std_out = sqrt(col_dev_means);
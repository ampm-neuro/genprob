function [pval, stats] = discrim_test_comb(all_paths, varargin)
%combines all session trl_mtx's and compares waits to rich and poor tones

% plot input
if nargin ==2
    plot_on = varargin{1};
else
    plot_on = 0;
end

% iterate through sessions combinging trl_mtx
trl_mtx_all = [];
for isesh = 1:size(all_paths,1)
    load(all_paths{isesh}, 'trl_mtx')
    trl_mtx_all = [trl_mtx_all; trl_mtx];
end

% compute wait times
[waits, tones] = wait_times_prep(trl_mtx_all, 2); 

% ttest
[~, pval, ~, stats] = ttest2(waits(tones==min(tones)), waits(tones~=min(tones)));

% plot
if plot_on > 0
    figure; 
    wait_times_plot(trl_mtx_all,2); 
end
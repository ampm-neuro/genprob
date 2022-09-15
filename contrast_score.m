function [num_pix] = contrast_score(act_vect)
% computes the minimum number of contiguous pixels holding 1/10 the value of vector


% minimum sum threshold
divisor = 10;
min_sum = sum(act_vect)/divisor;


%% any bins
%
% sort
sort_vect = sort(act_vect,'descend');

% cumsum
cum_sum = cumsum(sort_vect);

% find bin that passes threshold
num_pix_any = find(cum_sum>=min_sum, 1, 'first');
%}

%% contiguous bins
num_pix = [];
for i_bin_ct = num_pix_any:length(act_vect)
   for lo_bound = 1:length(act_vect)-(i_bin_ct-1)
       hi_bound = lo_bound + (i_bin_ct-1);
       if sum(act_vect(lo_bound:hi_bound))>=min_sum
           num_pix = i_bin_ct;
           return
       end
   end
end

    




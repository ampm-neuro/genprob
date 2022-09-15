function [wdw_mean, wdw_std, wdw_se] = nansmooth_ampm(vect_in, window_size)

wdw_mean = nan(length(vect_in), 1);
wdw_std = nan(size(wdw_mean));
wdw_se = nan(size(wdw_mean));

win_half = floor(window_size/2);
for islide = 1:length(vect_in)
    
    lo_bnd = islide-win_half; 
        if lo_bnd <1; lo_bnd=1; end
    hi_bnd = islide+win_half; 
        if hi_bnd >length(vect_in); hi_bnd=length(vect_in); end
        
        wdw_mean(islide) = nanmean(vect_in(lo_bnd:hi_bnd));
        wdw_std(islide) = nanstd(vect_in(lo_bnd:hi_bnd));
        wdw_se(islide) = nanstd(vect_in(lo_bnd:hi_bnd))/sqrt(sum(~isnan(vect_in(lo_bnd:hi_bnd))));
    
end
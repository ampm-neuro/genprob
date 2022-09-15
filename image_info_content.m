function info_content = image_info_content(fr_vect)
% computes info content of vector
    
% occupancy probabilities are flat
occ_p = repmat(1/length(fr_vect), 1, length(fr_vect));

%R is the overall mean firing rate
R = mean(fr_vect);

%calculate info content for each bin
info_contents = nan(size(fr_vect));
for ibin = 1:length(fr_vect)
    Pi = occ_p(ibin); % probability of occupancy of bin i
    Ri = fr_vect(ibin); % mean firing rate for bin i
    rr = Ri/R; % relative rate
    info_contents(ibin) = Pi*(rr)*log2(rr); % load info content for each bin
end



% sum for output
info_content = nansum(info_contents(:));
function pos_out = pdist_rwalk_2d(pdist, start_pos, num_steps)
% pos_out = pdist_rwalk_2d(pdist, start_pos, num_steps)
%
% outputs the sequential positions of a walk within the space
% defined by the probability distribution pdist (i.e. given an infinite
% number of steps, the histogram of output positions should perfectly 
% match pdist)
%
% at each step, the probability of each potential direction is computed as:
%
% the value for that adjacent bin over the value of all adjacent bins,
% where the value for an adjacent bin is computed as: 
% 
% the mean of the probability of the  adjacent bin (from pdist) and the 
% probability of the current bin.
%
% pdist is a vector of probabilities of occupying that each bin. 
% try normpdf(1:100, 50, 10)
%
% step size is the bin size of pdist
%
% the length of pos_out is num_steps + 1 for start position
% 

% pdist positions
pdist_pos = 1:length(pdist);

% preallocate positions
pos_out = nan(num_steps+1,1);
pos_out(1) = start_pos;

% iterate through each step
for istep = 1:length(pos_out)-1
    
    % select step direction
    p_left = round(sum(pdist([pos_out(istep) pos_out(istep)-1]))*100000);
    p_right = round(sum(pdist([pos_out(istep) pos_out(istep)+1]))*100000);
    step_direct = randsample([ones(p_left,1).*-1; ones(p_right,1)],1);
    
    % load
    pos_out(istep+1) = pos_out(istep) + step_direct;
    
end
    
    
end
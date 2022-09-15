function [exp_rwd_mtx, ideal_waits, waits] = expect_rwd(session_duration_s, bins, delay_distr, rwd_probs, varargin)
% plots a matrix of expected rewards for each mean wait duration (for rich
% and poor tones)
% 
% [exp_rwd_mtx, ideal_waits] = expect_rwd(60*60, 100, delay_distribution_gen4, [.9 .9 .1 .1], all_stage_learn_cell);

% load('delay_distribution_gen4')
% [all_stage_learn_mtx, all_stage_learn_cell] = ALL_stage_learn(1:11, 2, 1);


if nargin==5
    asl = varargin{1};
end

%fixed_delay
fixed_delay = 2; %s

%intertrial interval
iti = 1; %s

%waits
upper_wait_bound = 38; % max(delay_distr)
wait_times = linspace(0, upper_wait_bound, bins+1);
wait_times = mean([wait_times(1:end-1); wait_times(2:end)]);

%proportion of available rewards available at each wait time
prop_rwds = nan(length(wait_times),1);
for iwait = 1:length(wait_times)
    prop_rwds(iwait) = sum(delay_distr<wait_times(iwait))/length(delay_distr);
end

%dimensionality of output matrix
unq_prwd = unique(rwd_probs);
mtx_dim = length(unq_prwd);

%build matrix (hard coded for 2 dims)
exp_rwd_mtx = nan(length(wait_times), length(wait_times));
for iwait_dim1 = 1:length(wait_times)
    
    exp_rwd_dim1 = prop_rwds(iwait_dim1) * unq_prwd(1);
    
    for iwait_dim2 = 1:length(wait_times)
        
        
        exp_rwd_dim2 = prop_rwds(iwait_dim2) * unq_prwd(2);
        
        num_trials = session_duration_s / (mean(wait_times([iwait_dim1 iwait_dim2]))+fixed_delay+iti);
        
        mean_prop_rwds = mean([exp_rwd_dim1 exp_rwd_dim2]);
        
        exp_rwd_mtx(iwait_dim1, iwait_dim2) = num_trials * mean_prop_rwds;
        
    end
    
    
end

figure; imagesc(exp_rwd_mtx)
set(gca, 'ydir', 'normal');
set(gca,'TickLength',[0, 0]); box off;


%circle max
[i,j] = ind2sub(size(exp_rwd_mtx), find(exp_rwd_mtx==max(exp_rwd_mtx(:))));
hold on; plot(j,i, 'ro', 'markersize', 7)

ideal_rich_wait = wait_times(j) + fixed_delay;
ideal_poor_wait = wait_times(i) + fixed_delay;
ideal_waits = [ideal_rich_wait ideal_poor_wait];

xticklabels(wait_times(xticks) + fixed_delay);
yticklabels(wait_times(yticks) + fixed_delay);

 cbh = colorbar ; %Create Colorbar
 caxis([0 3*60])
 cbh.TickLabels = cbh.Ticks ./ 60 ;
 
 title('Rewards per minute')
 xlabel('Rich tone wait (s)')
 ylabel('Poor tone wait (s)')
 axis square
 hold on; plot([0 bins], [0 bins], 'k--')
 
 
%plot actual mouse sessions

asl_hold = cell(length(asl{1}),1);
%iterate through subjects
for i1 = 1:length(asl{1})
    keep_flag = 0;
    %iterate through stagesf
    for i2 = length(asl):-1:1
        if keep_flag ==1
            asl_hold{i1} = asl{i2}{i1,size(asl{i2},2)};
            break
        elseif ~isempty(asl{i2}{i1,1})
            keep_flag = 1;
        end
    end
end
asl = asl_hold;


rich_waits = []; poor_waits = []; 
for i = 1:length(asl)%[1 2 3 4 6 7 9:16]; 
    
    
    if isempty(asl{i})
        continue
    end
    
    rich_waits = [rich_waits; mean(asl{i}{1})]; 
    poor_waits = [poor_waits; mean(asl{i}{2})]; 
end
waits = [rich_waits, poor_waits];

for i = 1:size(waits,1)
    out_bins = waits_to_bins(wait_times, waits(i,:));
    hold on
    plot(out_bins(1), out_bins(2), 'ko', 'markersize', 7)
end
 
    function out_bins = waits_to_bins(bin_waits, raw_waits)
        bins_idx = 1:length(bin_waits);
        out_bins = nan(size(raw_waits));
        for iw = 1:length(raw_waits)
            out_bins(iw) = interp1(bin_waits, bins_idx, raw_waits(iw));
        end
    end



end
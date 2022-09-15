function [deviations] = p_dist_error_sim(num_sesh, num_round, shufs)
%calculate deviations from p_dist and cumulative reward history over n sessions

deviations = nan(shufs, num_sesh);

for ishuf = 1:shufs
    
    %simulate cumulative reward history
    [~, isesh_cum, p_dist] = p_dist_rwd_sim(num_sesh, num_round, 0);
    
    %compute deviations
    for icum = 1:size(isesh_cum,1)
    
        deviations(ishuf, icum) = sum(abs((p_dist./100)-isesh_cum(icum,:)));
    
    end
    
end

figure; hold on

plot(mean(deviations), 'k-')
plot(mean(deviations)+std(deviations), 'color', [1 1 1].*0.5)
plot(mean(deviations)-std(deviations), 'color', [1 1 1].*0.5)

%chance binary
plot(xlim, [1 1].*sum(abs(repmat(0.5,size(p_dist)) - p_dist./100)), 'k--')

set(gca,'TickLength',[0, 0]);
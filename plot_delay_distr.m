function plot_delay_distr(trl_mtx, fixed_delay)
% fixed delay is medass_cell{1}(6)/100

% REWARDS AVAILABLE
%plot delay durations (red are missed rewards)
bin_size = 2.0;
bin_bounds = 0:bin_size:60;
xtick_loc = linspace(min(bin_bounds), max(bin_bounds), length(bin_bounds))+diff(bin_bounds(1:2))/2; 
xtick_loc(end)=[];

obtained_rwd_delays = histcounts(trl_mtx(~isnan(trl_mtx(:,11)), 4), bin_bounds);
abandonded_rwd_delays = histcounts(trl_mtx(isnan(trl_mtx(:,11)) & trl_mtx(:,3)==1, 4), bin_bounds);

obtained_rwd_percent = sum(~isnan(trl_mtx(:,11)))/sum(trl_mtx(:,3)==1)





rwd_delay_mtx = [obtained_rwd_delays; abandonded_rwd_delays]';
rwd_delay_mtx = rwd_delay_mtx;
rwd_delay_mtx = rwd_delay_mtx./sum(rwd_delay_mtx(:)); %proportion
figure; bar(xtick_loc, rwd_delay_mtx, 'stacked')
set(gca,'TickLength',[0, 0]);
box off
xlim([0 60]); ylim([0 .25])
hold on; plot([1 1].*(fixed_delay), [0 1], 'k--')
legend({'Obtained', 'Abandoned'})
ylabel('Proportion of trials')
xlabel('Wait time (s)')

% REWARDS UNAVAILABLE
%plot delay durations
delay_durations = histcounts(trl_mtx(trl_mtx(:,3)==0,12)+ fixed_delay, bin_bounds) ;
delay_durations = delay_durations./sum(delay_durations(:));
figure; bar(xtick_loc, delay_durations)
set(gca,'TickLength',[0, 0]);
box off
xlim([0 60]); ylim([0 .25])
hold on; plot([1 1].*(fixed_delay), [0 1], 'k--')
ylabel('Proportion of trials')
xlabel('Wait time (s)')


%probes rich v poor
load('ALL_frequencies', 'all_unq_frq', 'rwd_idx')



obtained_rwd_percent_rich = sum(~isnan(trl_mtx(:,11)) & ismember(trl_mtx(:,2), all_unq_frq(rwd_idx)))/sum(trl_mtx(:,3)==1  & ismember(trl_mtx(:,2), all_unq_frq(rwd_idx)))
obtained_rwd_percent_poor = sum(~isnan(trl_mtx(:,11)) & ~ismember(trl_mtx(:,2), all_unq_frq(rwd_idx)))/sum(trl_mtx(:,3)==1  & ~ismember(trl_mtx(:,2), all_unq_frq(rwd_idx)))

rwd_percent_rich = sum(~isnan(trl_mtx(:,11)) & ismember(trl_mtx(:,2), all_unq_frq(rwd_idx)))/sum(ismember(trl_mtx(:,2), all_unq_frq(rwd_idx)))
rwd_percent_poor = sum(~isnan(trl_mtx(:,11)) & ~ismember(trl_mtx(:,2), all_unq_frq(rwd_idx)))/sum(~ismember(trl_mtx(:,2), all_unq_frq(rwd_idx)))


figure;hold on
    %poor
    delay_durations_poor = histcounts(trl_mtx(trl_mtx(:,3)==0 & ~ismember(trl_mtx(:,2), all_unq_frq(rwd_idx)),12)+ fixed_delay, bin_bounds) ;
    
    mean(trl_mtx(trl_mtx(:,3)==0 & ~ismember(trl_mtx(:,2), all_unq_frq(rwd_idx)),12)+ fixed_delay)
    
    delay_durations_poor = delay_durations_poor./sum(delay_durations_poor(:));
    bar(xtick_loc, delay_durations_poor)
    %rich
    delay_durations_rich = histcounts(trl_mtx(trl_mtx(:,3)==0 & ismember(trl_mtx(:,2), all_unq_frq(rwd_idx)),12)+ fixed_delay, bin_bounds) ;
    
    mean(trl_mtx(trl_mtx(:,3)==0 & ismember(trl_mtx(:,2), all_unq_frq(rwd_idx)),12)+ fixed_delay)
    delay_durations_rich = delay_durations_rich./sum(delay_durations_rich(:));
    bar(xtick_loc, delay_durations_rich)
plot(xtick_loc,smooth(delay_durations_poor))
plot(xtick_loc,smooth(delay_durations_rich))
set(gca,'TickLength',[0, 0]);
box off
xlim([0 60]); ylim([0 .25])
hold on; plot([1 1].*(fixed_delay), [0 1], 'k--')
ylabel('Proportion of trials')
xlabel('Wait time (s)')
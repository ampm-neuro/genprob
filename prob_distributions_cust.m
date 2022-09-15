h1 = figure; hold on
h2 = figure; hold on

binsize = .001;
num_delays = 500;
xponent = 2.00;
max_time = 60;
%{
fixed_delay = 0.25;
x_start = 0.25;
x = x_start:0.01:(x_start+max_time-fixed_delay); y=1./x.^xponent; y=y./sum(y);
g = 1000000; pdist = randpdf(y, x, [1 g])-x_start; 
pdist = sort(pdist); pdist = pdist(1:g/num_delays:end);
pdist = round(pdist.*100)./100;
figure(h1); histogram(pdist+fixed_delay, 0:binsize:max_time, 'normalization', 'probability')
figure(h2); plot(x - x_start + fixed_delay, y)
pdist_gen1 = pdist;
csvwrite('ddist_gen1.csv', pdist_gen1)

fixed_delay = 0.50;
x_start = 1.0;
x = x_start:0.01:(x_start+max_time-fixed_delay); y=1./x.^xponent; y=y./sum(y);
g = 1000000; pdist = randpdf(y, x, [1 g])-x_start; 
pdist = sort(pdist); pdist = pdist(1:g/num_delays:end);
pdist = round(pdist.*100)./100;
figure(h1); histogram(pdist+fixed_delay, 0:binsize:max_time, 'normalization', 'probability')
figure(h2); plot(x - x_start + fixed_delay, y)
pdist_gen2 = pdist;
csvwrite('ddist_gen2.csv', pdist_gen2)

fixed_delay = 1.00;
x_start = 2;
x = x_start:0.01:(x_start+max_time-fixed_delay); y=1./x.^xponent; y=y./sum(y);
g = 1000000; pdist = randpdf(y, x, [1 g])-x_start; 
pdist = sort(pdist); pdist = pdist(1:g/num_delays:end);
pdist = round(pdist.*100)./100;
figure(h1); histogram(pdist+fixed_delay, 0:binsize:max_time, 'normalization', 'probability')
figure(h2); plot(x - x_start + fixed_delay, y)
pdist_gen3 = pdist;
csvwrite('ddist_gen3.csv', pdist_gen3)
%}
%{
fixed_delay = 1.50;
x_start = 4;
x = x_start:1:(x_start+max_time-fixed_delay); y=1./x.^xponent; y=y./sum(y);
g = 1000000; pdist = randpdf(y, x, [1 g])-x_start; 
pdist = sort(pdist); pdist = pdist(1:g/num_delays:end);
pdist = round(pdist.*100)./100; 
figure(h1); histogram(pdist+fixed_delay, 0:binsize:max_time, 'normalization', 'probability')
figure(h2); plot(x - x_start + fixed_delay, y)
pdist_gen4 = pdist;
csvwrite('ddist_gen4.csv', pdist_gen4)
%}

%
fixed_delay = 2.00;
x_start = 10;
x = x_start:0.01:(x_start+max_time-fixed_delay); y=1./x.^xponent; y=y./sum(y);
g = 1000000; pdist = randpdf(y, x, [1 g])-x_start; 
pdist = sort(pdist); pdist = pdist(1:g/num_delays:end);
pdist = round(pdist.*100)./100; 
figure(h1); histogram(pdist+fixed_delay, 0:binsize:max_time, 'normalization', 'probability')
figure(h2); plot(x - x_start + fixed_delay, y)
pdist_gen5 = pdist;
csvwrite('ddist_gen5.csv', pdist_gen5)
%}


%delay index vector
del_idx = 0:num_delays-1;
csvwrite('delay_index_vector.csv', del_idx);
%}
%

%}
figure(h1);
set(gca,'TickLength',[0, 0]); ylim([0 1])

figure(h2);
set(gca,'TickLength',[0, 0]); axis auto
%close


% rich and poor
%{
figure; hold on
all_trl_mtx = ALL_sessions(1:10);
pdist_gen_cts = histcounts(all_trl_mtx(:,4)+fixed_delay, 0:binsize:max_time);
pdist_gen_cts_norm = pdist_gen_cts./sum(pdist_gen_cts(:));
pdist_rich = pdist_gen_cts_norm.*0.90;
pdist_poor = pdist_gen_cts_norm.*0.10;

pdist_revcum_rich = fliplr(cumsum(fliplr(pdist_rich)));
pdist_revcum_poor = fliplr(cumsum(fliplr(pdist_poor)));

bin_times = 0:binsize:max_time;
bin_times = mean([bin_times(1:end-1); bin_times(2:end)]);

plot(bin_times, fliplr(cumsum(fliplr(pdist_gen_cts_norm))), 'k-');
h1 = plot(bin_times, pdist_revcum_rich);
h2 = plot(bin_times, pdist_revcum_poor);

pdist_revcum = fliplr(cumsum(fliplr(pdist_gen_cts_norm.*0.5)));
pdist_cum = cumsum(pdist_gen_cts_norm);
%pdist_cum = 1-
h3 = plot(bin_times, pdist_revcum);

plot(bin_times, pdist_cum, 'k--')
plot(bin_times, pdist_cum.*0.90, 'k--', 'color', h1.Color)
plot(bin_times, pdist_cum.*0.50, 'k--', 'color', h3.Color)
plot(bin_times, pdist_cum.*0.42, 'k--', 'color', h3.Color)
plot(bin_times, pdist_cum.*0.10, 'k--', 'color', h2.Color)
%plot(pdist_cum.*.5, '--', 'color', h3.Color)
%plot(pdist_cum.*.9, '--', 'color', h1.Color)
%plot(pdist_cum.*.1, '--', 'color', h2.Color)


ylim([0 1])

%hold on; plot((12.9153).*[1 1], ylim, 'k-') %mean wait
hold on; plot((12.02).*[1 1], ylim, 'k-', 'color', h2.Color)
hold on; plot((20.667).*[1 1], ylim, 'k-', 'color', h1.Color)
hold on; plot(mean([12.02 20.667]).*[1 1], ylim, 'k-')

hold on; plot(xlim, [1 1].*0.50, 'k--')
hold on; plot(xlim, [1 1].*0.4194, 'k--')
hold on; plot(xlim, [1 1].*0.8526, 'k--', 'color', h1.Color)
hold on; plot(xlim, [1 1].*0.7279, 'k--', 'color', h2.Color)
hold on; plot(xlim, [1 1].*0.7644, 'k--', 'color', h1.Color)
hold on; plot(xlim, [1 1].*0.0730, 'k--', 'color', h2.Color)
set(gca,'TickLength',[0, 0]);
legend({'Rich p(rwd)', 'Poor p(rwd)'}, 'location', 'northeast')
%}

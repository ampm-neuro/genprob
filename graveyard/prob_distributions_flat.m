h1 = figure; hold on
h2 = figure; hold on

xponent = 2.00;
%
fixed_delay = 0.25;
x_start = 2;
x = x_start:0.01:(x_start+30-fixed_delay); y=1./x.^xponent; y=y./sum(y);
g = 1000000; pdist = randpdf(y, x, [1 g])-x_start; 
pdist = sort(pdist); pdist = pdist(1:g/100:end);
pdist = round(pdist.*100)./100;
figure(h1); histogram(pdist+fixed_delay, 0:1:30, 'normalization', 'probability')
figure(h2); plot(x - x_start + fixed_delay, y)
pdist_gen1 = pdist;
csvwrite('ddist_gen1.csv', pdist_gen1)
%}
fixed_delay = 0.50;
x_start = 5;
x = x_start:0.01:(x_start+30-fixed_delay); y=1./x.^xponent; y=y./sum(y);
g = 1000000; pdist = randpdf(y, x, [1 g])-x_start; 
pdist = sort(pdist); pdist = pdist(1:g/100:end);
pdist = round(pdist.*100)./100; 
figure(h1); histogram(pdist+fixed_delay, 0:1:30, 'normalization', 'probability')
figure(h2); plot(x - x_start + fixed_delay, y)
pdist_gen2 = pdist;
csvwrite('ddist_gen2.csv', pdist_gen2)
%
fixed_delay = 1.00;
x_start = 10;
x = x_start:0.01:(x_start+30-fixed_delay); y=1./x.^xponent; y=y./sum(y);
g = 1000000; pdist = randpdf(y, x, [1 g])-x_start; 
pdist = sort(pdist); pdist = pdist(1:g/100:end);
pdist = round(pdist.*100)./100; 
figure(h1); histogram(pdist+fixed_delay, 0:1:30, 'normalization', 'probability')
figure(h2); plot(x - x_start + fixed_delay, y)
pdist_gen3 = pdist;
csvwrite('ddist_gen3.csv', pdist_gen3)
%}
%
fixed_delay = 2.00;
x_start = 20;
x = x_start:0.01:(x_start+30-fixed_delay); y=1./x.^xponent; y=y./sum(y);
g = 1000000; pdist = randpdf(y, x, [1 g])-x_start; 
pdist = sort(pdist); pdist = pdist(1:g/100:end);
pdist = round(pdist.*100)./100; 
figure(h1); histogram(pdist+fixed_delay, 0:1:30, 'normalization', 'probability')
figure(h2); plot(x - x_start + fixed_delay, y)
pdist_gen4 = pdist;
csvwrite('ddist_gen4.csv', pdist_gen4)
%}
figure(h1);
set(gca,'TickLength',[0, 0]); ylim([0 1])

figure(h2);
set(gca,'TickLength',[0, 0]); axis auto
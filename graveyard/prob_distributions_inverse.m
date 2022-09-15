h1 = figure; hold on
h2 = figure; hold on

xponent = 2.00;
%
x_start = 0.25;
x = x_start:0.01:30; y=1./x.^xponent; y=y./sum(y);
g = 1000000; pdist = randpdf(y, x, [1 g]); 
pdist = sort(pdist); pdist = pdist(1:g/100:end);
pdist = round(pdist.*100)./100;
figure(h1); histogram(pdist, 0:1:30, 'normalization', 'probability')
figure(h2); plot(x, y)
pdist_gen1 = pdist-x_start;
csvwrite('pdist_gen1.csv', pdist_gen1)


x_start = 0.50;
x = x_start:0.01:30; y=1./x.^xponent; y=y./sum(y);
g = 1000000; pdist = randpdf(y, x, [1 g]); 
pdist = sort(pdist); pdist = pdist(1:g/100:end);
pdist = round(pdist.*100)./100; 
figure(h1); histogram(pdist, 0:1:30, 'normalization', 'probability')
figure(h2); plot(x, y)
pdist_gen2 = pdist-x_start;
csvwrite('pdist_gen2.csv', pdist_gen2)

x_start = 1.00;
x = x_start:0.01:30; y=1./x.^xponent; y=y./sum(y);
g = 1000000; pdist = randpdf(y, x, [1 g]); 
pdist = sort(pdist); pdist = pdist(1:g/100:end);
pdist = round(pdist.*100)./100; 
figure(h1); histogram(pdist, 0:1:30, 'normalization', 'probability')
figure(h2); plot(x, y)
pdist_gen3 = pdist-x_start;
csvwrite('pdist_gen3.csv', pdist_gen3)
%}
x_start = 2.00;
x = x_start:0.01:30; y=1./x.^xponent; y=y./sum(y);
g = 1000000; pdist = randpdf(y, x, [1 g]); 
pdist = sort(pdist); pdist = pdist(1:g/100:end);
pdist = round(pdist.*100)./100; 
figure(h1); histogram(pdist, 0:1:30, 'normalization', 'probability')
figure(h2); plot(x, y)
pdist_gen4 = pdist-x_start;
csvwrite('pdist_gen4.csv', pdist_gen4)

figure(h1);
set(gca,'TickLength',[0, 0]); ylim([0 1])

figure(h2);
set(gca,'TickLength',[0, 0]); ylim([0 1])
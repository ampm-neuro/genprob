
% randomish walk
pos_out = pdist_rwalk_2d(normpdf(1:100, 50, 10), 40, 100000);

%10s
figure; plot(pos_out(1:10), 1:10, '-o')
hold on; plot(mean(pos_out(1:10)).*[1 1], ylim, 'r-') %10s mean

%100s
figure; plot(pos_out(1:100), 1:100, '-o')
hold on; plot(mean(pos_out(1:10)).*[1 1], ylim, 'r-') %10s mean

%10000s
figure; plot(pos_out(1:10000), 1:10000, '-o')
hold on; plot(mean(pos_out(1:10)).*[1 1], ylim, 'r-') %10s mean
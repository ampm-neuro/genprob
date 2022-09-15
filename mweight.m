function new_amp = mweight(freqs)
% modify amplitude to produce equal perceived loudness in mice

% frequencies and approximate db thresholds for mice

% compute distance from intersection of mouse and 0
% add that value to volume measure in box
% may need to subtract from some

a = [5000 70;...
    6000 52;...
    7000 42;...
    10000 40;...
    15000 43;
    20000 46;...
    30000 56; ...
    40000 75];

figure; plot(a(:,1), a(:,2), 'o-');
set(gca, 'XScale', 'log')


new_amp = interp1(a(:,1), a(:,2), freqs, 'pchip');

hold on; plot(freqs, new_amp, 'o')


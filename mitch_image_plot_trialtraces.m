function mitch_image_plot_trialtraces(traces, neurons)

% normalize all traces to [0 1]
traces = traces(neurons,:);
traces = traces - min(traces, [], 2);
traces = traces ./ max(traces, [], 2);

% figure
figure; hold on; 

% iterate through all neurons
for ic = 1:length(neurons)
 
        % 1 cell per row
        yaxis_row = length(neurons) - ic + 0.5;

        % plot traces
        plot(traces(ic, :).*1.4 + yaxis_row);
end
  
% aesthetics
set(gcf,'Position', [550 55 1010 1283])
set(gca,'TickLength',[0, 0]); box off;
ylim([0.5 length(neurons)+0.5])
ytl_hold = length(neurons):-1:1;
yticklabels(ytl_hold(yticks))
ylabel('Neuron')
xlabel('Time (s)')
xlim([0 inf])
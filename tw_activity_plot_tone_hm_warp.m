function [all_matrices, xval, act_mtx] = tw_activity_plot_tone_hm_warp(trl_mtx, trl_idx, medass_cell, frame_times, traces, neurons, trials, align_event, tw_bounds, tses)
% runs tw_activity and produces one plot for each neuron comparing activity
% for each unique tone
%
% align events:
% NP onset = 1
% NP offset = 2
% Tone on = 3
% Head entry = 4
% Reward delivery = 5
% Reward receipt = 6
% Trial end = 7
%
% tw_bounds is a two item vector. the first item is the number of seconds
% before the align event the tw should start and the second item is the
% number of seconds after the align event the tw should end. eg [2 3] means
% the tw starts 2 seconds and ends 3s after the align event.
% compute all activity


%timewarp traces
[traces, frame_times, trl_idx, trl_mtx] = ...
    timewarp_traces(traces, frame_times, trl_idx, trl_mtx, tses);

%activity matrix
act_mtx = tw_activity_II(trl_mtx, trl_idx, frame_times, traces, neurons, trials, align_event, tw_bounds);

% constrain based on input
trl_mtx_local = trl_mtx(trials, :);

% unique tones
load('unqfrq41.mat', 'unqfrq41')
utones = unique(floor(trl_mtx(:,2)));

% soft preallocate
all_matrices = [];

% iterate through each neuron
for ic = 1:length(neurons)
    
    % iterate through tones
    tone_mtx = nan(length(unqfrq41),size(act_mtx,2));
    for itone = 1:length(utones)
        
        % compute mean and se
        trace_local = nanmean(act_mtx(ic, :, floor(trl_mtx_local(:,2))==utones(itone)), 3);
        %trace_local = norm_mtx(trace_local);
        tone_mtx(unqfrq41==utones(itone),:) = trace_local;

    end
    
    % plot
    figure; % needs to be off for hm_full
    hold on
    imagesc(tone_mtx)
    
    % aesthetics
    set(gca,'TickLength',[0, 0]); box off;
    ylim([0.5 length(unqfrq41)+0.5])
    xlim([0 inf])
    title(['Neuron number ' num2str(neurons(ic))])
    yticks(1:2:length(unqfrq41))
    yticklabels(unqfrq41(1:2:length(unqfrq41)))
    xticks([1 size(act_mtx,2)])
    xticklabels([-tw_bounds(1) tw_bounds(2)])
    xlabel('Time (s)')
    ylabel('Tone (Hz)')
    colorbar
    
    % red line
    xval = size(act_mtx,2).*(tw_bounds(1)/sum(tw_bounds)) - 0.5;
    plot(xval.*[1 1], ylim, 'r-')
    
    % load 3d matrix output
    all_matrices = cat(3, all_matrices, tone_mtx);
    
end

%% all matrix average
all_mtx = nanmean(all_matrices,3);

% plot
figure;
hold on
imagesc(all_mtx)

% aesthetics
set(gca,'TickLength',[0, 0]); box off;
ylim([0.5 length(unqfrq41)+0.5])
xlim([0 inf])
title(['Neuron number ' num2str(neurons(ic))])
yticks(1:2:length(unqfrq41))
yticklabels(unqfrq41(1:2:length(unqfrq41)))
xticks([1 size(act_mtx,2)])
xticklabels([-tw_bounds(1) tw_bounds(2)])
xlabel('Time (s)')
ylabel('Tone (Hz)')
colorbar

% red line
xval = size(act_mtx,2).*(tw_bounds(1)/sum(tw_bounds)) - 0.5;
plot(xval.*[1 1], ylim, 'r-')


%% line plots
%{
figure; hold on
colors = winter(size(all_mtx,2));
for i = 1:size(all_mtx,2)
    plot(all_mtx(:,i), '-.', 'color', colors(i,:))
end
%}
    
    
end
    
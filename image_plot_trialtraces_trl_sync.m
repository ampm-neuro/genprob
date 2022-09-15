function image_plot_trialtraces_trl_sync(trl_mtx, all_pkl_frames, trl_idx, traces, neurons, trials, varargin)
% plot each trial as a serperate figure; include all neurons sorted
% according to srt_idx (varargin)


min_size = 0.25; % s

% colors
%colors = [0 0 0; 186 32 13]./255;
%color_diff = (colors(2,:)-colors(1,:));

% remove abberations
for i = 1:size(traces)
    traces(i,:) = remove_abberations(traces(i,:), 2.5);
end

% sort neurons
if ~isempty(varargin)
   srt_idx = varargin{1};
else
    srt_idx = 1:length(neurons);
end
neurons = neurons(srt_idx);

% normalize all traces to [0 1]
traces = traces(neurons,:);
traces = traces - min(traces, [], 2);
traces = traces ./ max(traces, [], 2);

% iterate through all trials
for itrl = 1:length(trials)
    
    
    % figure
    figure; hold on; 
    
    % current items
    current_trial = trials(itrl);
    
    % time
    trl_time = all_pkl_frames(trl_idx==current_trial);
    min_time = min(trl_time);
    x_time = trl_time - min_time; % set start to 0s
    
    % id synronous activity
    trl_traces = traces(:, trl_idx==current_trial);
    trl_traces_sync = trl_traces-mean(trl_traces(:));
    
    % id higher than average population activity
    mean_all_traces = mean(trl_traces_sync);
    ste_lo = mean_all_traces - std(trl_traces_sync)./sqrt(size(trl_traces_sync,1));
    sync_idx = false(1, size(trl_traces_sync,2));
    for itw = 1:length(trl_time)
        current_time = trl_time(itw);
        tw_lo = current_time-min_size/2;
        tw_hi = current_time+min_size/2;
        time_idx = trl_time>=tw_lo & trl_time<=tw_hi;

        if all(ste_lo(time_idx)>0)
            sync_idx(time_idx)=true;
        end
    end
    
    % id higher than average individual activity
    %cuttoff_each_trace = max(trl_traces,[],2)./2;
    hi_cuttoff_each_trace = median(trl_traces,2)*1.25;
    lo_cuttoff_each_trace = median(trl_traces,2);
    indv_hi = false(size(trl_traces));
    for ic = 1:size(trl_traces,1)
        for itw = 1:length(trl_time)
            current_time = trl_time(itw);
            tw_lo = current_time-min_size/2;
            tw_hi = current_time+min_size/2;
            time_idx = trl_time>=tw_lo & trl_time<=tw_hi;
            
            if all(trl_traces(ic, time_idx)>hi_cuttoff_each_trace(ic)) && any(sync_idx(time_idx))
                indv_hi(ic, time_idx)=trl_traces(ic, time_idx)>lo_cuttoff_each_trace(ic);
            end
        end
    end

    
    % dual sync idx
    
    dual_sync_idx = indv_hi;% & repmat(sync_idx, size(indv_hi,1), 1);
    
    %figure; imagesc(repmat(sync_idx, size(indv_hi,1), 1))
    
    % sync plot prep
    ytrace_highlight = nan(size(trl_traces));
    ytrace_highlight(dual_sync_idx) = trl_traces(dual_sync_idx);
    
    
    % iterate through all neurons
    for ic = 1:length(neurons)
    
        % 1 cell per row
        yaxis_row = length(neurons) - ic + 0.5;
        
        % neuron activity
        y_trace = trl_traces(ic, :);
        
        % plot traces
        plot(x_time, y_trace.*1.5 + yaxis_row, 'color', .3.*[1 1 1]);

        % color plot traces
        %{
        z = zeros(size(x_time));
        col = y_trace;
        surface([x_time;x_time], [y_trace.*1.5 + yaxis_row; y_trace.*1.5 + yaxis_row], [z;z], [col;col],...
            'facecol','no','edgecol','interp','linew',1);
        %}
        
        % highlight syncronous traces
        plot(x_time, ytrace_highlight(ic,:).*1.5 + yaxis_row, 'b-', 'linewidth', 2);
        
      
        
        % plot NP exit
        np_offset = sum(trl_mtx(current_trial, [1 9])) - min_time;
        plot(np_offset.*[1 1], [0 1] + yaxis_row, 'r-', 'linewidth', 2)
        
        % plot head entry
        fd_onset = sum(trl_mtx(current_trial, [1 10])) - min_time;
        plot(fd_onset.*[1 1], [0 1] + yaxis_row, 'r-', 'linewidth', 2)
        
        % plot tone on
        tone_on_time = sum(trl_mtx(current_trial, [1 7])) - min_time;
        plot(tone_on_time.*[1 1], [0 1] + yaxis_row, 'r-', 'linewidth', 2)
        
        % plot start of random delay
        fd_offset = sum(trl_mtx(current_trial, [1 10])) - min_time + 2;
        plot(fd_offset.*[1 1], [0 1] + yaxis_row, 'r-', 'linewidth', 2)
     
        % plot reward time
        rwd_delivery_time = sum(trl_mtx(current_trial, [1 11])) - min_time;
        plot(rwd_delivery_time.*[1 1], [0 1] + yaxis_row, 'r-', 'linewidth', 2)
        
        % plot reward receipt time
        %{
        if ~isnan(sum(trl_mtx(current_trial, [1 11])))
            first_lick_time = min(medass_cell{15}(medass_cell{15}>sum(trl_mtx(current_trial, [1 11])))) - min_time;
            plot(first_lick_time.*[1 1], [0 1] + yaxis_row, '-', 'color', colors(itrl,:), 'linewidth', 2)
        end
        %}

    end
    
    % aesthetics
    set(gca,'TickLength',[0, 0]); box off;
    ylim([0.5 length(neurons)+0.5])
    ytl_hold = length(neurons):-1:1;
    yticklabels(ytl_hold(yticks))
    ylabel('Neuron')
    xlabel('Time (s)')
    xlim([0 max(x_time)])
    title(['Trial ' num2str(current_trial)])
    
end
  



% plot mean activity across all neurons
%
figure; hold on
trl_traces = traces(:, trl_idx==current_trial);
trl_traces = trl_traces-mean(trl_traces(:));

mean_trace = mean(trl_traces);
ste_lo = mean_trace - std(trl_traces)./sqrt(size(trl_traces,1));
ste_hi = mean_trace + std(trl_traces)./sqrt(size(trl_traces,1));

plot(x_time, mean_trace, 'k-', 'linewidth', 2)
plot(x_time, ste_lo, 'k-', 'linewidth', 1)
plot(x_time, ste_hi, 'k-', 'linewidth', 1)

plot(np_offset.*[1 1], [-1 1], 'r-', 'linewidth', 2)
plot(fd_onset.*[1 1], [-1 1], 'r-', 'linewidth', 2)
plot(tone_on_time.*[1 1], [-1 1], 'r-', 'linewidth', 2)
plot(fd_offset.*[1 1], [-1 1], 'r-', 'linewidth', 2)
plot(rwd_delivery_time.*[1 1], [-1 1], 'r-', 'linewidth', 2)

xlim([0 max(x_time)])
ylim([min(mean_trace)-.05 max(mean_trace)+.05])

% id synronous activity
trl_traces = traces(:, trl_idx==current_trial);
trl_traces = trl_traces-mean(trl_traces(:));
min_size = 1; %s
sync_idx = false(size(trl_traces));
for itw = 1:length(trl_time)
    current_time = trl_time(itw);
    tw_lo = current_time-min_size/2;
    tw_hi = current_time+min_size/2;
    time_idx = trl_time>=tw_lo & trl_time<=tw_hi;
    
    if all(ste_lo(time_idx)>0)
        sync_idx(time_idx)=true;
    end
end
hold on; plot(xlim, [0 0], 'k--')





% plot correlation between traces
%{
tw = 1; %s
trl_time = all_pkl_frames(trl_idx==current_trial);
trl_traces = traces(:, trl_idx==current_trial)';
all_cors_mean = nan(size(trl_time));
all_cors_std = nan(size(trl_time));
all_cors_se = nan(size(trl_time));
for itw = 1:length(trl_time)
    
    % time bounds
    current_time = trl_time(itw);
    tw_lo = current_time-tw/2;
    tw_hi = current_time+tw/2;
    time_idx = trl_time>=tw_lo & trl_time <=tw_hi;
    
    % all corrs
    rmtx = corrcoef(trl_traces(time_idx,:));
    
    % index for pairwise corrs
    if itw==1
        idx = false(size(rmtx));
        for irow = 1:size(rmtx,2)-1
            idx(irow,irow+1:end) = true;
        end
    end
    rmtx = rmtx(idx);
    
    % means etc
    all_cors_mean(itw) = mean(rmtx(:));
    all_cors_std(itw) = std(rmtx(:));
    
end

figure; hold on
plot(x_time, all_cors_mean, 'k-', 'linewidth', 2)
%plot(x_time, all_cors_mean+all_cors_std, 'color', 0.7.*[1 1 1])
%plot(x_time, all_cors_mean-all_cors_std, 'color', 0.7.*[1 1 1])
plot(x_time, all_cors_mean+all_cors_std./sqrt(length(rmtx(:))), 'color', 0.5.*[1 1 1])
plot(x_time, all_cors_mean-all_cors_std./sqrt(length(rmtx(:))), 'color', 0.5.*[1 1 1])

plot(np_offset.*[1 1], [-1 1], 'r-', 'linewidth', 2)
plot(fd_onset.*[1 1], [-1 1], 'r-', 'linewidth', 2)
plot(tone_on_time.*[1 1], [-1 1], 'r-', 'linewidth', 2)
plot(fd_offset.*[1 1], [-1 1], 'r-', 'linewidth', 2)
plot(rwd_delivery_time.*[1 1], [-1 1], 'r-', 'linewidth', 2)

xlim([0 inf])
ylim([-.2 .4])

%}






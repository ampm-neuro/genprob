function trace_footprint_activityPeak(footprint_mtx, trl_mtx, trl_idx, frame_times, traces, tses)


%identify peak of each cell
[~,srt_idx] = image_mean_activity_timewarp(trl_mtx, trl_idx, frame_times, traces, 1:size(traces,1), unique(trl_idx), tses);



% color according to relative position in session


% color according to ordinal position
colors = parula(length(srt_idx));
colors = colors(srt_idx,:);


% wrong number of cells
if size(footprint_mtx,1) ~= size(colors,1)
    error('mismatch')
end


% iterate through each cell
figure; hold on
for icell = 1:size(footprint_mtx,1)
    trace_footprint(squeeze(footprint_mtx(icell,:,:)), colors(icell,:));
end
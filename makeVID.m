% code for making gif files showing progression of single cell task
% responses over the trial

[a,srt_idx] = image_mean_activity_timewarp(...
    trl_mtx, trl_idx, frame_times, traces, 1:size(traces,1), unique(trl_idx));

desired_cells = srt_idx;

h = figure;
axis tight manual


% make new (to-be-written) video object
filename = 'trial_activity_prog';
new_vid_obj = VideoWriter(filename);
new_vid_obj.Quality = 98;
new_vid_obj.FrameRate = 22;

% open obj
open(new_vid_obj)

% details
%{
num_frames = nan(1, length(desired_cells));
fheights = nan(1, length(desired_cells));
fwidths = nan(1, length(desired_cells));
%}

for ifrm = 1:length(desired_cells)

    % make fig
    tw_activity_plot_trial_hm_warp(...
        trl_mtx, trl_idx, frame_times, traces, desired_cells(ifrm),...
        unique(trl_idx), [3.10 7.75], [0.2 1.1 1.5 2.0 3.5 2.0]);
    axis off; colorbar off; title([])
    set(gcf, 'Position', [1352 337 560 654])
    drawnow
    
    % Capture the plot as an image 
    frame = getframe(h); 
    %im = frame2im(frame); 
    %[imind,cm] = rgb2ind(im,256);
    
    % write video
    writeVideo(new_vid_obj, frame);
    
end

% close obj
close(new_vid_obj)
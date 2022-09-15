% code for making gif files showing progression of single cell task
% responses over the trial

[a,srt_idx] = image_mean_activity_timewarp(...
    trl_mtx, trl_idx, frame_times, traces, 1:size(traces,1), unique(trl_idx));

desired_cells = srt_idx;

h = figure;
axis tight manual
filename = 'trial_activity_prog.gif';

for ifrm = 1:length(desired_cells)

    % make fig
    tw_activity_plot_trial_hm_warp(...
        trl_mtx, trl_idx, frame_times, traces, desired_cells(ifrm),...
        unique(trl_idx), [3.5 8], [0.2 1.1 1.5 2.0 3.5 2.0]);
    axis off; colorbar off; title([])
    drawnow 
    
    % Capture the plot as an image 
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256);
    
    % initiate gif
    if ifrm == 1
        imwrite(imind,cm,filename,'gif','DelayTime',0.1, 'Loopcount',inf); 
    else
        imwrite(imind,cm,filename,'gif','DelayTime',0.1,'WriteMode','append'); 
    end
    
    
    
    
end
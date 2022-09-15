function simul_vid(vid_paths, varargin)
% outputs a single video that simultaneously plays each of the videos 
% listed in vid_paths
%
% vid_paths is a cell containing the paths to each to-be-simul-played
% video

% argin
if nargin == 2
    frame_times = varargin{1};
else
    frame_times = cell(length(vid_paths),1);
end


% number of videos
num_vids = length(vid_paths);

% make video objects for each video path
vid_objs = cell(1,num_vids);
for ivid = 1:num_vids
    vid_objs{ivid} = VideoReader(vid_paths{ivid});
end

% num frames
num_frames = nan(1, num_vids);
fheights = nan(1, num_vids);
fwidths = nan(1, num_vids);
for ivid = 1:num_vids
    num_frames(ivid) = vid_objs{ivid}.NumberOfFrames;
    fheights(ivid) = vid_objs{ivid}.Height;
    fwidths(ivid) = vid_objs{ivid}.Width;
   
    % frame times
    if isempty(frame_times{ivid})
        frame_times{ivid} = 1:num_frames(ivid);
    end
end


% frame time info
ft_min = min(cellfun(@min, frame_times)); % min frametime
ft_max = max(cellfun(@max, frame_times)); % max frametime
fr_max = min(cellfun(@(x) min(x(2:end)-x(1:end-1)), frame_times)); % min time between frames



% iterate through time, updating video frames when they occur
for itime = ft_min : fr_max : ft_max
    
    % iteratively update each subplot axis with current frame
    for ivid = 1:num_vids
        subplot(1, num_vids, ivid);

        % find current frames        
        ft_min = abs(frame_times{ivid} - itime);
        if min(ft_min)<.25
            iframe = find(ft_min==min(ft_min));
            vid_frm = vid_objs{ivid}.read(iframe);
        else
            vid_frm = zeros(fheights(ivid), fwidths(ivid), 3, 'uint8');
        end
        image(vid_frm);
        axis off
    end
    
    % update figure displays
    set(gcf,'Position', [336 791 1391 425])
    drawnow
    
end








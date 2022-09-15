function simul_vid_write(vid_paths, new_vid_filename, varargin)
% outputs a single video that simultaneously plays each of the videos 
% listed in vid_paths
%
% vid_paths is a cell containing the paths to each to-be-simul-played
% video
%
% for info on frame time input, see the function simul_vid_sync

% argin
if nargin == 3
    frame_times = varargin{1};
else
    frame_times = cell(length(vid_paths),1);
end


% number of videos
num_vids = length(vid_paths);

% make video objects for each video path
vid_objs = cell(1,num_vids);
im_idx = false(1,num_vids);
for ivid = 1:num_vids
    if contains(vid_paths{ivid}, '.tif')
        %vid_objs{ivid} = VideoReader(vid_paths{ivid});
        im_idx(ivid) = true;
    else
        vid_objs{ivid} = VideoReader(vid_paths{ivid});
    end
end

% make new (to-be-written) video object
new_vid_obj = VideoWriter(new_vid_filename);
new_vid_obj.Quality = 98;

% num frames
num_frames = nan(1, num_vids);
fheights = nan(1, num_vids);
fwidths = nan(1, num_vids);
for ivid = 1:num_vids
    
    if im_idx(ivid) == 1
        im_info = imfinfo(vid_paths{ivid});
        fheights_hold = [im_info.Height];
        fheights(ivid) = fheights_hold(1);
        fwidths_hold = [im_info.Width];
        fwidths(ivid) = fwidths_hold(1);
        num_frames(ivid) = length(fwidths_hold);
    else
        num_frames(ivid) = vid_objs{ivid}.NumberOfFrames;
        fheights(ivid) = vid_objs{ivid}.Height;
        fwidths(ivid) = vid_objs{ivid}.Width;
    end
    
    % frame times
    if isempty(frame_times{ivid})
        frame_times{ivid} = 1:num_frames(ivid);
    end

end


% frame time info
ft_min = min(cellfun(@min, frame_times)); % min frametime
ft_max = max(cellfun(@max, frame_times)); % max frametime
fr_max = min(cellfun(@(x) min(x(2:end)-x(1:end-1)), frame_times)); % min time between frames

% update to-be-written video object
new_vid_obj.FrameRate = round(1/fr_max);

%figure
h1 = figure;

% open obj
open(new_vid_obj)

% limit to frames surrounding mm vid


% iterate through time, updating video frames when they occur
fig_handles = cell(1,2);
ft_min = 30; ft_max = 61.5;
last_vids = cell(1,num_vids);
for itime = ft_min : fr_max : ft_max
    
    % iteratively update each subplot axis with current frame
    for ivid = 1:num_vids
        fig_handles{ivid} = subplot(1, num_vids, ivid);

        % find current frames   
        ft_min = abs(frame_times{ivid} - itime);
        if min(ft_min)<.25
            iframe = find(ft_min==min(ft_min));
            if im_idx(ivid) == 1
                vid_frm = imread(vid_paths{ivid},iframe);
                vid_frm = vid_frm(95:395,230:535,:);
                %vid_frm = imadjust(vid_frm);
            else
                vid_frm = vid_objs{ivid}.read(iframe);
            end
        else
            vid_frm = zeros(fheights(ivid), fwidths(ivid), 3, 'uint8');
        end
        
        
        % zoom in on illuminated portion
        if ivid == 2
            
        end
        
        % update figure displays
        set(gcf,'Position', [336 791 900 425])
        if num_vids==2 && ivid == 1
            pos1 = fig_handles{1}.Position;
            pos1_new = [0.03 pos1(2) pos1(3)+0.15 pos1(4)];
            set(fig_handles{1},'Position', pos1_new)
            image(vid_frm);
        elseif num_vids==2 && ivid == 2
            pos2 = fig_handles{2}.Position;
            pos2 = [sum(pos1_new([1 3])+0.01) pos2(2) pos2(3)+0.1 pos2(4)];
            set(fig_handles{2},'Position', pos2)
            imagesc(vid_frm); %axis square
            colormap gray
            caxis([40 260])
        end
        axis off

    end

    drawnow
    
    % capture as frame
    frame = getframe(h1); 
    
    % write video
    writeVideo(new_vid_obj, frame);
    
    
end

% close obj
close(new_vid_obj)








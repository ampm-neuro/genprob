function [short_vid_file, OnOff_times_new, OnOff_times_odl, lum_mean] = gen_vid_preproc(vid_fp)
% takes a raw mkv video file and removes frames when LED was off
% outputs new vid file and matrices of timesstamps (for new and old vid) 
% with times of the begining and end of each LED period

% preallocate luminence vector
videoObject = VideoReader(vid_fp);
numberOfFrames = ceil(videoObject.Duration * videoObject.FrameRate);
lum_mean = nan(numberOfFrames,1);

% iterate through frames
videoObject = VideoReader(vid_fp);

frame_num = 0;
while hasFrame(videoObject)
    frame_num = frame_num+1;
    
    % read frame
    img = readFrame(videoObject);
    
    % grayscale image
    img = rgb2gray(img);
    
    % mean luminance
    lum_mean(frame_num) = mean2(img);
end

figure; plot(lum_mean)

short_vid_file = [];
OnOff_times_new = [];
OnOff_times_odl = [];

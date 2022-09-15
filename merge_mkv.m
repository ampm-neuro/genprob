function merge_mkv(video_fp)
% merges all mkv files in the folderpath video_fp

% get all mkv paths
direct = dir(video_fp);
direct = struct2cell(direct);
all_paths = cellfun(@(x1, x2) [x1 '\' x2], direct(2,:), direct(1,:), 'UniformOutput', false)';
videoList = all_paths(contains(all_paths, '.mkv') & ~contains(all_paths, 'mkv.pkl'));
videoList = sort(videoList);

% create output in seperate folder (to avoid accidentally using it as input)
outputVideo = VideoWriter([video_fp, '\mergedVideo'], 'Grayscale AVI');

% if all clips are from the same source/have the same specifications
% just initialize with the settings of the first video in videoList
inputVideo_init = VideoReader(videoList{1}); % first video
outputVideo.FrameRate = inputVideo_init.FrameRate;

% open stream
open(outputVideo)

% iterate over all videos you want to merge (e.g. in videoList)
for i = 1:length(videoList)
    
    disp(['Working on video ' num2str(i) ' of ' num2str(length(videoList))])
    
    % select i-th clip (assumes they are in order in this list!)
    inputVideo = VideoReader(videoList{i});
    
    % stream your inputVideo into an outputVideo
    while hasFrame(inputVideo)
        frm = readFrame(inputVideo);
        frm = rgb2gray(frm);
        writeVideo(outputVideo, frm);
    end
    
end

% close after having iterated through all videos
close(outputVideo)

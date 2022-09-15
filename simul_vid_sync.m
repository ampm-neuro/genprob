function frame_times = simul_vid_sync(vid_paths, frame_times_v1)
% computes the frame times to syncronize playback via simul_vid function
% is hard coded to deal with a medass video and a minimic video with
% corresponding pkl_frame_times input


%vid_paths{1} = 'D:\658648m2\gen09\d02\Box 1_04Mar2020_11-46-01.wmv';
%vid_paths{2} = 'C:\Users\ampm1\Desktop\mmm\mergedVideo.avi';

% frame_times_v1 is the frame_time output from the pkl files


% preallocate
frame_times = cell(length(vid_paths),1);

% med ass video
v_ma = VideoReader(vid_paths{1});
v_ma_fr = (1 / v_ma.Framerate);
v_ma_num = v_ma.NumberofFrames;

%0 : v_ma_fr :(v_ma_num*v_ma_fr)

frame_times{1} = 0 : v_ma_fr :(v_ma_num*v_ma_fr);


% minimic video
frame_times{2} = frame_times_v1;

% play videos
%simul_vid(vid_paths, frame_times)

% write videos
simul_vid_write(vid_paths, [fileparts(vid_paths{2}) '\simulVid.avi'], frame_times)



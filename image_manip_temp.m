rotate_angle = 25; % counterclockwise degrees
xshift = 100; yshift = 20; % pixels
translate_shift = [xshift yshift];


% load image
%foldername = 'D:\Projects\InProgress\Mitch_LinTrack\control\152-2\2021-04-06';
%mov_file = get_file_paths_targeted(foldername, {'.avi'});
mov_file= 'C:\Users\ampm1\Documents\vid_test\2021-03-10-15-38-44.mkv';
%mov = VideoReader(mov_file{1});
mov = VideoReader(mov_file);
im1 = readFrame(mov);
%im1 = im1+100 ;


figure; hold on
h1 = imagesc(im1);

figure; hold on
size(im1)
size(repmat(reshape(rgb2lab([1 1 1]), 1, 1, 3), size(im1,1), size(im1,2), 1))
h1 = imagesc(double(im1).*repmat(reshape(rgb2lab([1 1 1]), 1, 1, 3), size(im1,1), size(im1,2), 1));
 
%set(h1, 'AlphaData', [.5]);

%im2 = im1;
%im2 = imrotate(im2, rotate_angle);
%im2 = imtranslate(im2, translate_shift,'OutputView','full');

%h2 = image(im2);
%set(h2, 'AlphaData', [.5]);

function [all_pkl_frames, trl_idx] = process_splitvids(video_fp)
% combines minimic videos and pkl files recorded during a single medass
% session
% input path to split video and pkl files, and path to raw medass file

% medass file is always first
medass_fp = dir(video_fp);
try
medass_fp = [video_fp '\' medass_fp(1).name];
catch
    dir(video_fp)
    [video_fp '\' medass_fp(1).name]
    
end

%% delete invalid video files
delete_tag = 0;
[mkv_paths, pkl_paths] = delete_invalid_mkv(video_fp, medass_fp); delete_tag=1;

% get all paths
%
direct = dir(video_fp);
direct = struct2cell(direct);
all_paths = cellfun(@(x1, x2) [x1 '\' x2], direct(2,:), direct(1,:), 'UniformOutput', false)';

% get all mkv paths
mkv_paths = all_paths(contains(all_paths, '.mkv') & ~contains(all_paths, 'mkv.pkl'));

% get all pkl paths
pkl_paths = all_paths(contains(all_paths, '.pkl') & ~contains(all_paths, '.csv'));
%}

%% merge remaining pkl files
all_pkl_frames = [];
trl_idx = [];
for ipklpath = 1:length(pkl_paths)

    % path components
    [filepath,name,ext] = fileparts(pkl_paths{ipklpath});
    pkl_folder = filepath;
    pkl_file = [name,ext];
    
    % convert to csv then to vector
    pkl_frames_local = loadpickle(pkl_folder, pkl_file);
    
    % merge pkl frames
    all_pkl_frames = [all_pkl_frames; pkl_frames_local];
    
    % add trial index
    trl_idx = [trl_idx; repmat(ipklpath, size(pkl_frames_local))];
end
%figure; plot(trl_idx, all_pkl_frames)


% delete csv files
direct = dir(video_fp);
direct = struct2cell(direct);
all_paths = cellfun(@(x1, x2) [x1 '\' x2], direct(2,:), direct(1,:), 'UniformOutput', false)';
csv_files_all = all_paths(contains(all_paths, '.csv'));
for itbd = csv_files_all'
    delete(itbd{:})
end

%% synch pkl timestamps with med associates file

%
% load medass
[medass_cell, trl_mtx] = load_medass(medass_fp);

% plot trials
%figure; plot_trials_trlmtx(trl_mtx)


%{
for itrl = unique(trl_idx)'

    
    if itrl > length(medass_cell{2})
        break
    end
    
    
    % first video on (in medass time) coincides first valid trial
    vid_start_time_ma = medass_cell{2}(itrl); % trial start

    % pkl start of first video time
    sv_pkl = min(all_pkl_frames(trl_idx==itrl));

    % correct pkl times
    all_pkl_frames(trl_idx==itrl) = all_pkl_frames(trl_idx==itrl) + (vid_start_time_ma - sv_pkl);

end
%}
%
vid_start_time_ma = medass_cell{2}(1); %trial start
sv_pkl = min(all_pkl_frames);
all_pkl_frames = all_pkl_frames + (vid_start_time_ma - sv_pkl);
%}




%%  plot
%
figure; hold on
ts_ma = sum(trl_mtx(:,[1 6]),2); % trial start medass (np onset)
te_ma = sum(trl_mtx(:,[1 12]),2); % trial end medass
ts_pkl = nan(length(unique(trl_idx)),1); % trial start pkl
te_pkl = nan(length(unique(trl_idx)),1); % trial end pkl


for itrl = 1:200 %medass_cell{1}(29) %unique(trl_idx)'     

    % plot medass TRIAL IS GREY
    try
        patch([ts_ma(itrl) te_ma(itrl) te_ma(itrl) ts_ma(itrl)], [0 0 .6 .6], 0.5.*[1 1 1], 'FaceAlpha',0.5)
        text(ts_ma(itrl),-.05, num2str(itrl))
    catch
    end
    % plot pkl VIDEO IS green if valid, red if invalid
    try
            ts_pkl(itrl) = min(all_pkl_frames(trl_idx==itrl));
            te_pkl(itrl) = max(all_pkl_frames(trl_idx==itrl));
            if delete_tag==0 && medass_cell{17}(itrl) == 0
                patch([ts_pkl(itrl) te_pkl(itrl) te_pkl(itrl) ts_pkl(itrl)], [.4 .4 1 1], 'r', 'FaceAlpha',0.5)
            else
                patch([ts_pkl(itrl) te_pkl(itrl) te_pkl(itrl) ts_pkl(itrl)], [.4 .4 1 1], 'g', 'FaceAlpha',0.5)
            end
            text(ts_pkl(itrl),1.05, num2str(itrl))
    catch
    end
end
ylim([-.1 1.1])
%}


%% save new combined csv file
save([video_fp '\pkl_frametimes.mat'], 'all_pkl_frames')
save([video_fp '\pkl_trl_idx.mat'], 'trl_idx')


%% delete original pkl files
%{
pkl_files_all = all_paths(contains(all_paths, '.pkl') & ~contains(all_paths, '.csv'));
for itbd = pkl_files_all'
    delete(itbd{:})
end
%}


%% merge remaining mkv files
% does not work with cnmfe for some reason. use python based solution
%merge_mkv(video_fp)


%% delete original mkv files
%{
mkv_files_all = all_paths(contains(all_paths, '.mkv') & ~contains(all_paths, 'mkv.pkl'));
for itbd = mkv_files_all'
    delete(itbd{:})
end
%}


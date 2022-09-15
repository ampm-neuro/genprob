function [mkv_paths, pkl_paths] = delete_invalid_mkv(video_fp, medass_fp)
% uses index in medass_cell 17 (Q) to erase video files from invalid trials
%(0s in the index)

% get valid trials index
[medass_cell, ~] = load_medass(medass_fp);

% get all paths
direct = dir(video_fp);
direct = struct2cell(direct);
all_paths = cellfun(@(x1, x2) [x1 '\' x2], direct(2,:), direct(1,:), 'UniformOutput', false)';

% get all mkv paths
mkv_files_all = all_paths(contains(all_paths, '.mkv') & ~contains(all_paths, 'mkv.pkl'));

% get all pkl paths
pkl_files_all = all_paths(contains(all_paths, '.pkl') & ~contains(all_paths, '.csv'));

%
disp('NumVideos NumTrialInitations NumTrialCompletions NumValidVideos')
[length(mkv_files_all) length(medass_cell{17}) medass_cell{1}(16) sum(medass_cell{17})]
%{
if length(mkv_files_all) ~= length(medass_cell{17})
    if length(mkv_files_all) ~= medass_cell{1}(16)
        error
    end
end
%}

% check if process was previously run
if length(mkv_files_all) ~= length(medass_cell{17})
    mkv_paths = mkv_files_all;
    pkl_paths = pkl_files_all;
    
    if length(mkv_files_all) == sum(medass_cell{17}) && exist([video_fp '\discard'], 'dir')
        disp('Incompatible medass index and video list. Invalid videos removed previously.')
        return
    elseif length(mkv_files_all) == length(medass_cell{17})-1 %extra unfinished trial
        disp('Last medass_cell{17} item ignored.')
        medass_cell{17} = medass_cell{17}(1:end-1);
    else
        disp('Incompatible medass index and video list.')
        video_fp
        error
        return
    end
end

% identify invalid and valid mkv file paths
to_be_deleted_paths_mkv = mkv_files_all(~logical(medass_cell{17}));
mkv_paths = mkv_files_all(logical(medass_cell{17}));

% identify invalid and valid pkl file paths
to_be_deleted_paths_pkl = pkl_files_all(~logical(medass_cell{17}));
pkl_paths = pkl_files_all(logical(medass_cell{17}));

% "discard" invalid files
%
for itbd = [to_be_deleted_paths_mkv; to_be_deleted_paths_pkl]'
    
    % create discard folder
    if ~exist([video_fp '\discard'], 'dir')
        mkdir([video_fp '\discard'])
    end
    
    % move file to discard folder
    movefile(itbd{:}, [video_fp '\discard'])
end
%}

%check if destination has a containing folder for old files.
%{
if ~exist([video_fp '/discard'], 'dir')
    mkdir([video_fp '/discard'])
end
%}

%{
allpaths = [mkv_paths; pkl_paths];
for itbm = allpaths
    [a,b,c] = fileparts(allpaths{itbm});
    movefile(allpaths, [a '\old\' b c])
end
%}



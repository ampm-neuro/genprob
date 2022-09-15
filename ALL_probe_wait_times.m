function [awt_mtx, unique_subjs, subj_path_cell] = ALL_probe_wait_times(training_group, varargin)
% raw wait times from every probe for every subject (problem,tone,subject)

% smooth style
% 0 = raw (with nans)
% 1 = interp and exterp
% 2 = interp, exterp, and smooth


%% input
if nargin > 1
    smooth_style = varargin{1};
else
    smooth_style = 0; %raw with nans
end


%% file prep
% paths
folder_in = ['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\' training_group];
all_probe_paths = get_file_paths_targeted(folder_in, 'probe');
%all_probe_paths = all_probe_paths(~contains(all_probe_paths, 'notone'));
all_probe_paths = all_probe_paths(~contains(all_probe_paths, 'quiet'));

% unique subjects
all_subj_ids = [];
for ipath = 1:size(all_probe_paths,1)
    all_subj_ids = [all_subj_ids; find_subj_id(all_probe_paths{ipath})];
end
unique_subjs = unique(all_subj_ids, 'rows');

% one cell per subject
subj_path_cell = cell(1, size(unique_subjs,1));
for isubj = 1:size(unique_subjs,1)
    subj_path_cell{isubj} = all_probe_paths(contains(all_probe_paths, unique_subjs(isubj,:)) & contains(all_probe_paths, 'pre'));
    subj_path_cell{isubj} = [subj_path_cell{isubj}; all_probe_paths(contains(all_probe_paths, unique_subjs(isubj,:)) & contains(all_probe_paths, 'post') & ~contains(all_probe_paths, 'notone'))];
    no_tone_paths = all_probe_paths(contains(all_probe_paths, unique_subjs(isubj,:)) & contains(all_probe_paths, 'notone'));
    
    if ~isempty(no_tone_paths)
        subj_path_cell{isubj} = [subj_path_cell{isubj}; no_tone_paths(size(no_tone_paths,1))];
    end
end


%% compute wait times

% preallocate output matrix (probe, tone, subject)
awt_mtx = nan(8,41,size(unique_subjs,1));

% iterate through subject cells
for isubj = 1:length(subj_path_cell)
    
    % iterate through subject paths
    for isubj_path = 1:size(subj_path_cell{isubj})
        
        % load session
        load(subj_path_cell{isubj}{isubj_path}, 'trl_mtx')
        
        % compute and load wait times
        mean_or_all = 1; % means
        data_smoothing_protocol = 0; % model fit
        awt_mtx(isubj_path, :, isubj) = wait_times_prep(trl_mtx, mean_or_all, data_smoothing_protocol);
        %awt_mtx(isubj_path, :, isubj) = log(wait_times_prep(trl_mtx, mean_or_all, data_smoothing_protocol));
    end 
end


%% smoothing, interpolation, extrapolation

% smooth each session
%
if smooth_style == 1 % interp and extrap
    awt_mtx_hold = nan(size(awt_mtx));
    for isubj = 1:length(subj_path_cell)
        for isubj_path = 1:size(subj_path_cell{isubj})
            if sum(~isnan(awt_mtx(isubj_path, :, isubj)))>0
                %awt_mtx_hold(isubj_path, :, isubj) = nansmooth_adaptive_ampm(awt_mtx(isubj_path, :, isubj), 3, 3);
            end
        end 
    end 
    awt_mtx(isnan(awt_mtx)) = awt_mtx_hold(isnan(awt_mtx));
elseif smooth_style == 2 % interp, extrap, AND smooth
    for isubj = 1:length(subj_path_cell)
        for isubj_path = 1:size(subj_path_cell{isubj})
            if sum(~isnan(awt_mtx(isubj_path, :, isubj)))>0
                awt_mtx(isubj_path, :, isubj) = nansmooth_adaptive_ampm(awt_mtx(isubj_path, :, isubj), 3, 3);
            end
        end
    end 
end
%}











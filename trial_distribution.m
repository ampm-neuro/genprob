function trl_starts = trial_distribution(varargin)
% plots distribution of trial start times over the course of the session

% umbrella path
folderpath = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\Gen_richards';

% constrain via string inputs
fpaths = get_file_paths_targeted(folderpath, varargin);

% load trial start times
trl_starts = [];
for ifile = 1:length(fpaths)
   
    % load session
    load(fpaths{ifile}, 'trl_mtx')
    
    % trial start times
    trl_starts = [trl_starts; trl_mtx(:,1)];
    
end

% set to minutes
trl_starts = trl_starts./60;

% distribution
trl_start_dist = histcounts(trl_starts, 0:1:60);


%figure; 
%plot(trl_start_dist)

trl_start_dist = trl_start_dist./sum(trl_start_dist);
plot(cumsum(trl_start_dist))

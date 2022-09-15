function [sigmas] = plot_probes_one_subj(datafolder, subject)
% overlays all probe files that contain the input string

% folderpath
folderpath = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\';
folderpath = [folderpath datafolder];
folderpath = [folderpath '\' subject];

%title prep
fpaths_title = get_file_paths_targeted(folderpath, 'newlearn0')
if isempty(fpaths_title)
    NL = [];
    NL_str = 'NL';
else
    NL = fpaths_title{1}(strfind(fpaths_title{1}, 'newlearn0') + length('newlearn0'));
    NL_str = ['NL' NL];
end

% filter with additional input string constraints
probe_str = 'probe';
fpaths = get_file_paths_targeted(folderpath, probe_str)
if isempty(fpaths)
    error('No matching files found')
else
    delay = fpaths{1}(strfind(fpaths{1}, probe_str) + length(probe_str) + 1 : strfind(fpaths{1}, probe_str) + length(probe_str)+2);
    delay_str = ['D' delay];
end

% colors
colors = distinguishable_colors(length(fpaths));

% legend prep
figure; hold on
legend_input = cell(length(fpaths),1);
for ifile = 1:length(fpaths)
   plot(6500,1,'-', 'color', colors(ifile,:), 'linewidth', 1.5)
   legend_input{ifile} = ['p' num2str(ifile)];
end

% preallocate output
gof_stats = nan(length(fpaths),1);
pvals = nan(length(fpaths),1);
sigmas = nan(length(fpaths),1);

% iteratively load and plot files
for ifile = 1:length(fpaths)
    
    %load
    load(fpaths{ifile}, 'trl_mtx', 'medass_cell');
    
    %zscore waits
    %trl_mtx(trl_mtx(:,3)==0,12) = zscore_mtx(trl_mtx(trl_mtx(:,3)==0,12));
    
    % plot waits
    plot_normal_fit_subj(trl_mtx, 1, colors(ifile,:));
    
    % compute NORMAL goodness of fit stats
    [coefEsts] = plot_normal_fit_subj(trl_mtx, 2, colors(ifile,:));
    coefEsts
    sigmas(ifile,1) = coefEsts(4);


end

%update legend with stats
for ili = 1:length(legend_input)
    legend_input{ili} = [legend_input{ili} 'sigma=' num2str(round((sigmas(ili)*100))/100)];
end

% NL loc
if ~isempty(NL)
    NL_tone_lines(str2num(NL));
end

% legend
legend(legend_input, 'location', 'northeastoutside')

% title
title([NL_str '; ' delay_str])


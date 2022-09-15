function [wait_times_all] = plot_allprobes_ctx(folder_path, probe_number)
% plot indicated preprobes; Black for first preprobe and grey for second
% preprobe (black is generally in the new context, grey in the old)

%folder_path = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_6ctx';

probe_number
% identify all relevant files
all_first_probe_paths = [];
all_last_probe_paths = [];
for iprobe = 1:length(probe_number)
    if probe_number < 7
        all_probe_paths = [get_file_paths_targeted(folder_path, ['preprobe_0' num2str(probe_number(iprobe))])];
    elseif probe_number >=7
        all_probe_paths = [get_file_paths_targeted(folder_path, ['postprobe_0' num2str(probe_number(iprobe)-6)])];
    end

    % unique subjects
    all_subjs = [];
    for ipath = 1:size(all_probe_paths,1)
        all_subjs = [all_subjs; find_subj_id(all_probe_paths{ipath})];
    end
    all_subjs = unique(all_subjs, 'rows');

    % organize them by subject (1 cell per subj)
    subj_cell = cell(1, size(all_subjs,1));
    for isubj = 1:size(all_subjs,1)
        subj_cell{isubj} =  all_probe_paths(contains(all_probe_paths, all_subjs(isubj,:)));
    end

    % all first and last probe paths
    for isubj = 1:size(all_subjs,1)
        all_first_probe_paths = [all_first_probe_paths; subj_cell{isubj}(1)];
        all_last_probe_paths = [all_last_probe_paths; subj_cell{isubj}(end)];
    end

end

% plot
colors = [.3;.7].*[1 1 1;1 1 1];
figure; hold on; 
 

all_first_probe_paths

wait_times_all = cell(1,2);
[~,~,wait_times_all{1}] = plot_allprobes_colorin(all_first_probe_paths, colors(1,:));
[~,~,wait_times_all{2}] = plot_allprobes_colorin(all_last_probe_paths, colors(2,:));

%{
hold on; 
if probe_number < 7
    plot_allprobes_colorin(get_file_paths_targeted('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar\', {['preprobe_0' num2str(probe_number)]}), [.5 1 0])
elseif probe_number >=7
    plot_allprobes_colorin(get_file_paths_targeted('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar\', {['postprobe_0' num2str(probe_number)-6]}), [.5 1 0])
end

hold on; 
if probe_number < 7
    plot_allprobes_colorin(get_file_paths_targeted('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_hivar\', {['preprobe_0' num2str(probe_number)]}), [0 1 1])
elseif probe_number >=7
    plot_allprobes_colorin(get_file_paths_targeted('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_hivar\', {['postprobe_0' num2str(probe_number)-6]}), [0 1 1])
end
%}

try
rich_bounds_prob('mevar', probe_number-1);
catch
end
    
function all_probe_smooth_meta(region_string, probe_nums)
% runs plot_allprobes_opto_fit for every probe and then plots
% coefficients


% colors
colors = [0 0 1; .5 .5 .5]; %distinguishable_colors(2);

% preallocation
coefficient_cell = cell(1,length(probe_nums));
subject_cell = cell(1,length(probe_nums));

%% all probes smoothed
probe_ct = 0;
for iprobe = probe_nums
    probe_ct = probe_ct+1;
    if iprobe <=6
        path_str_1 = 'preprobe_opto';
        path_str_2 = ['_0' num2str(iprobe) '_'];
    elseif iprobe > 6 && iprobe < 9
        path_str_1 = 'postprobe_opto';
        path_str_2 = ['_0' num2str(iprobe-6) '_'];
    elseif iprobe == 9
        path_str_1 = 'postprobe_notone_opto';
        path_str_2 = '_02_';
    end
    paths = get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_optoExp_', region_string], path_str_1, path_str_2);
    
    %paths = paths(3)
    
    % subject cell
    for ipath = 1:size(paths,1)
        subject_cell{probe_ct} = [subject_cell{probe_ct}; find_subj_id(paths{ipath})];
    end
    
    % coefficient cell
    figure; hold on; 
    [~, coefficient_cell{probe_ct}] = plot_allprobes_opto_fit(paths);
    ylim([0 50])
    title(['Probe ' num2str(iprobe)])
    if iprobe>1
        rich_bounds_prob('mevar', iprobe-1, 1); 
    end
end


%% subjects
all_unique_subjects = [];
for iprobe = 1:length(subject_cell)
    all_unique_subjects = [all_unique_subjects; subject_cell{iprobe}];
end
all_unique_subjects = unique(all_unique_subjects, 'rows');


%% coefficients
for iprobe = 1:length(subject_cell)
    nan_hold = nan(size(all_unique_subjects, 1), size(coefficient_cell{iprobe}, 2), 2);

    if size(subject_cell{iprobe},1)>size(unique(subject_cell{iprobe},'rows'),1)
        [~,idx] = unique(subject_cell{iprobe},'rows', 'last');
        coefficient_cell{iprobe} = coefficient_cell{iprobe}(idx, :, :);
    end

    nan_hold(ismember(all_unique_subjects, subject_cell{iprobe}, 'rows'), :, :) = coefficient_cell{iprobe};
    coefficient_cell{iprobe} = nan_hold;
    
end


%% plot each coefficient across all probes
plot_cell = cell(size(coefficient_cell));

titles = cell(1,4);
titles{1} = 'Intercept';
titles{2} = 'Inflection freq';
titles{3} = 'Logistic';
titles{4} = 'Normal';

% one plot per coefficient
for icoef = [1 3 4]
    
    % figure
    figure; hold on
    
    % one errorbar for each laser condition
    for laser = 1:2
        
        if icoef == 1
            % intercept is summed with logistic to give poor-side min
            for iprobe = 1:length(coefficient_cell)
                plot_cell{iprobe} = coefficient_cell{iprobe}(:,1, laser) + coefficient_cell{iprobe}(:,2, laser);
            end
        elseif icoef == 2
            % flip based on rich side
            for iprobe = 1:length(coefficient_cell)
                if rem(iprobe,2)==0
                    plot_cell{iprobe} = -coefficient_cell{iprobe}(:,icoef, laser);
                else
                    plot_cell{iprobe} = coefficient_cell{iprobe}(:,icoef, laser);
                end
            end
        else
            for iprobe = 1:length(coefficient_cell)
                plot_cell{iprobe} = coefficient_cell{iprobe}(:,icoef, laser);
            end
        end
        
        % plot
        errorbar_plot(plot_cell, 1, [], colors(laser,:))
        
    end
    
    % aesthetics
    title(titles{icoef})
    
end


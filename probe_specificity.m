function probe_specificity(training_group, probe_num, varargin)
% probe_specificity(training_group, probe_num, problem_num, number_of_adj_tones)
% measure how much the probe reflects the rich and poor tones of the last
% training problem
%
% can input training problem number to use those tones, otherwise defaults 
% to the training problem before selected probe
%
% can input number of adjacent tones to compare to target, otherwise 
% defaults to one on each side
%
% probe must be 2:8, as probe 1 (pre-probe) has no associated training 
% problem

%% check inputs
if ~isempty(varargin)
    problem_num = varargin{1}; 
    if isempty(problem_num)
        problem_num = probe_num-1;
    end
    if length(varargin) == 2
        ajct_tns = varargin{2};
    end
else
    problem_num = probe_num-1;
    ajct_tns = 1;
end


%% compare rich and poor tones relative to adjacent tones
 
% get relevant probe sessoions
training_folder = ['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_' training_group];
if probe_num <= 6
    path_str = ['preprobe_0' num2str(probe_num)];
elseif probe_num >= 7
    path_str = ['postprobe_0' num2str(probe_num-6)];
end
all_paths = get_file_paths_targeted(training_folder, path_str);

% preallocate waits
rich_waits = nan(size(all_paths,1), 2*ajct_tns + 1);
poor_waits = nan(size(rich_waits));

% get rich and poor tone freqs
if contains(training_folder, 'mevar')
    rp_tones = rich_bounds_prob('mevar', 0);
elseif contains(training_folder, 'hivar')
    rp_tones = rich_bounds_prob('hivar', 0);
end
if probe_num <= 7
    rp_tones = rp_tones(problem_num, :);
elseif probe_num >= 8
    rp_tones = rp_tones(end, :);    
end

% for each probe session
for ipath = 1:size(all_paths,1)

    % load session
    load(all_paths{ipath})
    
    % compute waits
    [~, frequencies_all] = wait_times_prep(trl_mtx, 1, 1); %probe
    wait_times = nan(size(frequencies_all));
    [wait_times_hold, frequencies] = wait_times_prep(trl_mtx, 1, 0); %probe
    wait_times(ismember(frequencies_all, frequencies)) = wait_times_hold;
    frequencies = frequencies_all;
    
    if length(frequencies) ~= 41
        error
    end
    
    % rich and poor tone numbers
    rich_tone_num = find(frequencies==rp_tones(1));
    poor_tone_num = find(frequencies==rp_tones(2));
    
    % find upper and lower bounds
    rlo = rich_tone_num-ajct_tns; 
        ajct_tns_rlo_load = ajct_tns;
        if rlo < 1; ajct_tns_rlo_load = ajct_tns - (1-rlo); rlo = 1; end
        rlo_load = ceil(size(rich_waits,2)/2) - ajct_tns_rlo_load;
    rhi = rich_tone_num+ajct_tns; 
        ajct_tns_rhi_load = ajct_tns;
        if rhi > 41; ajct_tns_rhi_load = ajct_tns - (rhi - 41); rhi = 41; end
        rhi_load = ceil(size(rich_waits,2)/2) + ajct_tns_rhi_load;
        
    plo = poor_tone_num-ajct_tns; 
        ajct_tns_plo_load = ajct_tns;
        if plo < 1; ajct_tns_plo_load = ajct_tns - (1-plo); plo = 1; end
        plo_load = ceil(size(poor_waits,2)/2) - ajct_tns_plo_load;
    phi = poor_tone_num+ajct_tns; 
        ajct_tns_phi_load = ajct_tns;
        if phi > 41; ajct_tns_phi_load = ajct_tns - (phi - 41); phi = 41; end
        phi_load = ceil(size(poor_waits,2)/2) + ajct_tns_phi_load;
    
    % load wait times
    rich_waits(ipath, rlo_load : rhi_load) = wait_times(rlo : rhi);
    poor_waits(ipath, plo_load : phi_load) = wait_times(plo : phi);

end

% zscore waits
rich_waits = zscore_mtx(rich_waits')';
poor_waits = zscore_mtx(poor_waits')';

%% plot

% prepare errorbar plots
rw_cell = [];
for i = 1:size(rich_waits,2)
    rw_cell = [rw_cell {rich_waits(:,i)}];
end

pw_cell = [];
for i = 1:size(poor_waits,2)
    pw_cell = [pw_cell {poor_waits(:,i)}];
end

% plot waits by freq
figure;

subplot(2,2,1)
hold on
errorbar_plot(rw_cell, 1); 
title rich
plot(ceil(size(poor_waits,2)/2) .* [1 1], ylim, 'r-')

subplot(2,2,2)
hold on; 
errorbar_plot(pw_cell, 1); 
title poor
plot(ceil(size(poor_waits,2)/2) .* [1 1], ylim, 'r-')

ajct_tns_test = 1;
mid_point = ceil(length(rw_cell)/2);
test_lo = mid_point - ajct_tns_test;
test_hi = mid_point + ajct_tns_test;

subplot(2,2,3)
hold on
errorbar_plot([{nanmean([rw_cell{test_lo} rw_cell{test_hi}],2)} rw_cell(ceil(length(rw_cell)/2))], 1); 
title rich
[~, rich_pval] = ttest(mean([rw_cell{test_lo} rw_cell{test_hi}],2), rw_cell{ceil(length(rw_cell)/2)})

subplot(2,2,4)
hold on; 
errorbar_plot([{nanmean([pw_cell{test_lo} pw_cell{test_hi}],2)} pw_cell(ceil(length(pw_cell)/2))], 1); 
title poor
[~, poor_pval] = ttest(mean([pw_cell{test_lo} pw_cell{test_hi}],2), pw_cell{ceil(length(pw_cell)/2)})




%%
% fit line to group data, compare rich and poor tone relative to line



% compare rich and poor tones relative to mean wait times
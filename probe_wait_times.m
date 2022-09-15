function [mean_waits, subj_ids, all_waits] = probe_wait_times(folder_in, probes_in, tones_in, varargin)
% plot the average wait time over all tones_in for each probes_in

% folder to get probe session files from
%folder_in = 'mevar_meta';
dfp = ['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\' folder_in];

% preallocate
mean_waits = cell(1,length(probes_in));
all_waits = cell(1,length(probes_in));
subj_ids = cell(1,length(probes_in));

% dprime?
if ~isempty(varargin)
    raw_or_dprime = varargin{1};
else
    raw_or_dprime = 1; % difference
end

% iterate through probes
for iprobe = probes_in

    if iprobe <= 6

        % probe string
        ps = ['_0' num2str(iprobe) '_'];
    
        ppaths = get_file_paths_targeted(dfp, 'preprobe', ps);
        ppaths = ppaths(~contains(ppaths, 'notone'));
        ppaths = ppaths(~contains(ppaths, 'quiet'));
        
        % compute probe wait times
        [~, ~, all_wait_times] = noplot_allprobes(ppaths, 0);
    
    elseif iprobe > 6

        % probe string
        ps = ['_0' num2str(iprobe-6) '_'];
    
        ppaths = get_file_paths_targeted(dfp, 'postprobe', ps);
        ppaths = ppaths(~contains(ppaths, 'quiet'));
        ppaths = ppaths(~contains(ppaths, 'notone'));
        
        % compute probe wait times
        [~, ~, all_wait_times] = noplot_allprobes(ppaths, 0);
    end
    
    if isempty(all_wait_times)
        continue
    end
    
    % unique subjs
    usubjs = [];
    for ipath = 1:size(ppaths,1)
        usubjs = [usubjs; find_subj_id(ppaths{ipath})];
    end
    subj_ids{iprobe} = usubjs;
    
    % first subj ids
    if iprobe == 1
        first_subj_ids = usubjs;
    end
    shared_subjects = ismember(first_subj_ids, usubjs, 'rows');
    
    % insert nans for missing subjects
    all_wait_times_hold = all_wait_times;
    all_wait_times = nan(size(shared_subjects,1), size(all_wait_times_hold,2));
    all_wait_times(shared_subjects, :) = all_wait_times_hold;    
    
    
    % load average wait times
    %mean_waits{iprobe} = nanmean(all_wait_times(:,tones_in),2);

    all_waits{iprobe} = all_wait_times(:,tones_in);
    
    % how to compute wait times?
    mean_waits{iprobe} = nan(size(shared_subjects));
    if raw_or_dprime == 0
    	mean_waits{iprobe}(shared_subjects) = nanmean(all_wait_times(shared_subjects,tones_in),2);
    elseif raw_or_dprime == 1
        mean_waits{iprobe}(shared_subjects) = nanmean(all_wait_times(shared_subjects,tones_in),2) ./ nanmean(all_wait_times(shared_subjects, setdiff(1:41, tones_in)),2);
    elseif raw_or_dprime == 2
        %for irw = 1:size(all_wait_times,1)
        for irw = find(shared_subjects==1)'
            %mean_waits{iprobe}(irw,1) = zdiff(all_wait_times(irw,tones_in), all_wait_times(irw, setdiff(1:41, tones_in)));
            awt_toi = all_wait_times(irw,tones_in); awt_toi = awt_toi(:);
            awt_ntoi = all_wait_times(irw,setdiff(1:41, tones_in)); awt_ntoi = awt_ntoi(:);
            mean_waits{iprobe}(irw) = computeCohen_d(awt_toi, awt_ntoi);
        end
    elseif raw_or_dprime == 3
        mean_waits{iprobe}(shared_subjects) = nanmean(all_wait_times(shared_subjects,tones_in),2) - nanmean(all_wait_times(shared_subjects, setdiff(1:41, tones_in)),2);
    end
    
    %all_waits{iprobe} = all_wait_times(:,tones_in) ./ all_wait_times;
    
    
end

% plot
line_plot = nan(length(probes_in));
hold on

for pm = length(mean_waits):-1:1
    if isempty(mean_waits{pm})
        try
            mean_waits(pm) = [];
            subj_ids(pm) = [];
            all_waits{iprobe}(pm) = [];
        catch
        end
    end
end
if contains(folder_in, 'hivar')
    errorbar_plot(mean_waits, 1, [], .8.*[1 1 1])
else
    errorbar_plot(mean_waits, 1, [], .3.*[1 1 1])
end
%{
ct = 0;
for iprobe = probes_in
    ct = ct+1;
    line_plot(ct) = mean(probe_means{iprobe});
    errorbar(iprobe, mean(probe_means{iprobe}), std(probe_means{iprobe})./sqrt(length(probe_means{iprobe})), 'k')
end

% line
plot(probes_in, line_plot, 'k')
%}

% aesthetics
%ylim([0 30])
xlim([0 iprobe+1])
set(gca,'TickLength',[0, 0]); box off;





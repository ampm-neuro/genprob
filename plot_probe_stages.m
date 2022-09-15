function [all_wait_times_out, all_subj_out] = plot_probe_stages(probe_stages, folder_string)
% plots selected probe stages

%figure; 
hold on; 

dfp = ['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\' folder_string];

all_wait_times_out = [];
all_subj_out = [];

% iterate through probes
%iprobe_ct = 0;
for iprobe = probe_stages
    %iprobe_ct = iprobe_ct +1;
    

    if iprobe <= 6

        % probe string
        ps = ['_0' num2str(iprobe) '_'];
    
        % compute probe wait times
        ppaths = get_file_paths_targeted(dfp, 'preprobe', ps);
        ppaths = ppaths(~contains(ppaths, {'notone', 'quiet'}));
        [~, ~, all_wait_times] = plot_allprobes(ppaths);
        xlim([4500 38500])
        
        if iprobe==2
            rich_bounds_prob('mevar', 1);
        end
    
    elseif iprobe > 6 && iprobe <= 8

        % probe string
        ps = ['_0' num2str(iprobe-6) '_'];
    
        % compute probe wait times
        ppaths = get_file_paths_targeted(dfp, 'postprobe', ps);
        ppaths = ppaths(~contains(ppaths, {'notone', 'quiet'}));
        [~, ~, all_wait_times] = plot_allprobes(ppaths);
        xlim([4500 38500])
        
        if iprobe==7
            rich_bounds_prob('mevar', 6);
        end
        
    elseif iprobe ==9
    
        % compute probe wait times
        ppaths = get_file_paths_targeted(dfp, 'postprobe_notone');
        %ppaths = ppaths(end);
        [~, ~, all_wait_times] = plot_allprobes(ppaths);
        xlim([4500 38500])
        
    end

    all_wait_times_out = [all_wait_times_out; all_wait_times];
    
    % subjects
    for ipath = 1:size(ppaths,1)
        all_subj_out = [all_subj_out; repmat(find_subj_id(ppaths{ipath}), size(all_wait_times))];
    end
    
    mawt = nanmean(all_wait_times);
    mawt = 13.0;
    hold on; plot(xlim, mawt.*[1 1], 'k--')

    
    % normal fit
    %{
    [all_trl_mtx] = ALL_trl_mtx(ppaths);
    plot_normal_fit_subj(all_trl_mtx, 2);
    %}
    
end



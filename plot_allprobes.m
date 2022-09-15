function [all_fixed_effects, all_beta_pvals, wait_times_all] = plot_allprobes(file_paths)
% plot all probes
% varargin is colors

green_color = [21 154 72]./255; blue_color =  [43 73 159]./255;
light_green_color = [162 201 129]./255; light_blue_color = [102 164 213]./255;
if contains(file_paths{1}, 'mevar')
   primary_color = green_color;
   secondary_color = light_green_color;
else
    primary_color = blue_color;
   secondary_color = light_blue_color;

end

hold on

if isempty(file_paths)
    all_fixed_effects = [];
    all_beta_pvals = [];
    wait_times_all = [];
    return
end

% identify probe numbers
pp_nums = nan(size(file_paths,1),1);
for ipath = 1:size(file_paths,1)
   
    u_idx = strfind(file_paths{ipath}, '_');
    probe_num_string = file_paths{ipath}(u_idx(end-1)+1 : u_idx(end)-1);
    
    pp_nums(ipath) = str2double(probe_num_string);
    
end



% preallocate
all_fixed_effects = [];
all_beta_pvals = [];


for inum = 1:max(pp_nums)
    
    ppnum_fps = file_paths(pp_nums==inum);

    %plot means
    %
    if ~isempty(ppnum_fps)
        
        %get wait times from every session
        wait_times_all = [];
        frequencies_all = [];
        for ippnum = 1:size(ppnum_fps,1)
            load(ppnum_fps{ippnum}, 'trl_mtx')
            [wait_times_local, frequencies_local] = wait_times_prep(trl_mtx, 1, 0);
            if size(wait_times_local,2)>size(wait_times_local,1)
                wait_times_local = wait_times_local';
            end
            
            hold on; plot(frequencies_local, wait_times_local, '-', 'color', secondary_color)
            
            %wait_times_local = zscore_mtx(wait_times_local);
            %wait_times_local = smooth(wait_times_local, 3);
            wait_times_all = [wait_times_all; wait_times_local];
            frequencies_all = [frequencies_all; frequencies_local];
        end
        
        % means and ses
        unq_frq = unique(frequencies_all);
        mean_wait_times = nan(length(unq_frq),1);
        se_wait_times = nan(size(mean_wait_times));
        for ifrq = 1:length(unq_frq)
            mean_wait_times(ifrq) = nanmean(wait_times_all(frequencies_all==unq_frq(ifrq)));
            se_wait_times(ifrq) = nanstd(wait_times_all(frequencies_all==unq_frq(ifrq)))/sqrt(sum(frequencies_all==unq_frq(ifrq)));
        end
        
        % plot
       
        
        %color_hold = rand(1,3);
        errorbar(unq_frq, mean_wait_times, se_wait_times, 'linewidth', 1.5, 'color', primary_color)
        %errorbar(unq_frq, mean_wait_times, se_wait_times, 'linewidth', 1.5)
        set(gca,'TickLength',[0, 0]);
        xlabel('Tone Frequency (Hz)')
        ylabel('Wait Durations (s)')
        ylim_hold = ylim; ylim([0 ylim_hold(2)]);
        set(gca, 'XScale', 'log')
        title('Waits')
        xlim([4500 42000])
        xticks([5000 8500 14000 23000 35000])
        
        
        %[all_trl_mtx] = ALL_trl_mtx(ppnum_fps);
        %wait_times_plot(all_trl_mtx,1,colors(inum,:))
        
        
        
    else
        %plot([1 1+realmin],[1 1], 'color', colors(inum,:))
    end
    %}
    
    
end
ylim([0 60])


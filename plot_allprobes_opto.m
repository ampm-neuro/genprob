function [all_fixed_effects, all_beta_pvals] = plot_allprobes_opto(file_paths)
% plot all preprobes


hold on
file_paths
%get_file_paths_targeted('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\repeat_probes', {'preprobe', 'consistent'})

% opto tones
load('unqfrq41.mat', 'unqfrq41')
tones_optoON = unqfrq41(1:2:41);
tones_optoOFF = unqfrq41(2:2:40);

% identify probe numbers
pp_nums = nan(size(file_paths,1),1);
for ipath = 1:size(file_paths,1)
   
    u_idx = strfind(file_paths{ipath}, '_');
    probe_num_string = file_paths{ipath}(u_idx(end-1)+1 : u_idx(end)-1);
    
    pp_nums(ipath) = str2double(probe_num_string);
    
end

hold on
colors = [0 0 1; .5 .5 .5]; %distinguishable_colors(2);

% legend
legend_trick(colors, '-')

% preallocate
all_fixed_effects = [];
all_beta_pvals = [];

% iterate through unique probe numbers
for inum = 1:max(pp_nums)
    
    ppnum_fps = file_paths(pp_nums==inum); 
    %}
    %ppnum_fps = file_paths;

    %plot means
    %
    if ~isempty(ppnum_fps)
        
        %get wait times from every session
        wait_times_all_optoON = nan(1, length(unqfrq41));
        frequencies_all_optoON = nan(1, length(unqfrq41));
        wait_times_all_optoOFF = nan(1, length(unqfrq41));
        frequencies_all_optoOFF = nan(1, length(unqfrq41));
        for ippnum = 1:size(ppnum_fps,1)
            load(ppnum_fps{ippnum}, 'trl_mtx')
            
            figure; wait_times_plot_opto(trl_mtx,3); ylim([0 50])
            title(ppnum_fps{ippnum})
            
            % opto ON
            [wait_times_local_optoON, frequencies_local_optoON] = wait_times_prep(trl_mtx(ismember(floor(trl_mtx(:,2)), tones_optoON),:), 1, 0);
            
            wt_hold = nan(1, length(unqfrq41));
            wt_hold(ismember(unqfrq41, frequencies_local_optoON)) = wait_times_local_optoON;
            fl_hold = nan(1, length(unqfrq41));
            fl_hold(ismember(unqfrq41, frequencies_local_optoON)) = frequencies_local_optoON;
            
            wait_times_all_optoON = [wait_times_all_optoON; wt_hold];
            frequencies_all_optoON = [frequencies_all_optoON; fl_hold];
            
            % opto OFF
            [wait_times_local_optoOFF, frequencies_local_optoOFF] = wait_times_prep(trl_mtx(ismember(floor(trl_mtx(:,2)), tones_optoOFF),:), 1, 0);
            
            wt_hold = nan(1, length(unqfrq41));
            wt_hold(ismember(unqfrq41, frequencies_local_optoOFF)) = wait_times_local_optoOFF;
            fl_hold = nan(1, length(unqfrq41));
            fl_hold(ismember(unqfrq41, frequencies_local_optoOFF)) = frequencies_local_optoOFF;
            
            wait_times_all_optoOFF = [wait_times_all_optoOFF; wt_hold];
            frequencies_all_optoOFF = [frequencies_all_optoOFF; fl_hold];
            
        end
        
        % means and ses
        mean_wait_times_optoON = nanmean(wait_times_all_optoON);
        se_wait_times_optoON = nanstd(wait_times_all_optoON)./sqrt(sum(~isnan(wait_times_all_optoON),1));
        mean_wait_times_optoOFF = nanmean(wait_times_all_optoOFF);
        se_wait_times_optoOFF = nanstd(wait_times_all_optoOFF)./sqrt(sum(~isnan(wait_times_all_optoOFF),1));
        
        % plot
        figure; hold on
        errorbar(unqfrq41(~isnan(mean_wait_times_optoON)), mean_wait_times_optoON(~isnan(mean_wait_times_optoON)), se_wait_times_optoON(~isnan(mean_wait_times_optoON)), 'linewidth', 1.5, 'color', colors(1,:))
        errorbar(unqfrq41(~isnan(mean_wait_times_optoOFF)), mean_wait_times_optoOFF(~isnan(mean_wait_times_optoOFF)), se_wait_times_optoOFF(~isnan(mean_wait_times_optoOFF)), 'linewidth', 1.5, 'color', colors(2,:))
        
        
        
        
        
        %aesthetics
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
    %else
    %    plot([1 1+realmin],[1 1], 'color', colors(inum,:))
    end
    %}
end

ylim([0 60])
rich_bounds_2
legend({'LaserON', 'LaserOFF'}, 'location', 'northeast')
%rich_bounds
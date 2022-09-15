function [mean_wait_times_optoON, std_wait_times_optoON, se_wait_times_optoON, mean_wait_times_optoOFF, std_wait_times_optoOFF, se_wait_times_optoOFF] = plot_allprobes_opto_smooth(file_paths, varargin)
% plot all preprobes

%get_file_paths_targeted('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\repeat_probes', {'preprobe', 'consistent'})

if nargin==2
    plot_on = varargin{1};
else
    plot_on = 1;
end

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

if plot_on >= 1
    hold on
    colors = [0 0 1; .5 .5 .5]; %distinguishable_colors(2);
    legend_trick(colors, '-')
end

% iterate through unique probe numbers
for inum = max(pp_nums)
    
   % ppnum_fps = file_paths(pp_nums==inum); 
    %}
    ppnum_fps = file_paths;

    %plot means
    %
    if ~isempty(ppnum_fps)
        
        %get wait times from every session
        wait_times_all_optoON = [];
        frequencies_all_optoON = [];
        wait_times_all_optoOFF = [];
        frequencies_all_optoOFF = [];
        for ippnum = 1:size(ppnum_fps,1)
            load(ppnum_fps{ippnum}, 'trl_mtx')
            
            %figure; wait_times_plot_opto(trl_mtx,3); ylim([0 50])
            %title(ppnum_fps{ippnum})
            
            % opto ON
            [wait_times_local_optoON, frequencies_local_optoON] = wait_times_prep(trl_mtx(ismember(floor(trl_mtx(:,2)), tones_optoON),:), 1, -1);
            
            wt_hold = nan(1, length(unqfrq41));
            wt_hold(ismember(unqfrq41, frequencies_local_optoON)) = wait_times_local_optoON;
            fl_hold = nan(1, length(unqfrq41));
            fl_hold(ismember(unqfrq41, frequencies_local_optoON)) = frequencies_local_optoON;
            
            wait_times_all_optoON = [wait_times_all_optoON; wt_hold];
            frequencies_all_optoON = [frequencies_all_optoON; fl_hold];
            
            % opto OFF
            [wait_times_local_optoOFF, frequencies_local_optoOFF] = wait_times_prep(trl_mtx(ismember(floor(trl_mtx(:,2)), tones_optoOFF),:), 1, -1);
            
            wt_hold = nan(1, length(unqfrq41));
            wt_hold(ismember(unqfrq41, frequencies_local_optoOFF)) = wait_times_local_optoOFF;
            fl_hold = nan(1, length(unqfrq41));
            fl_hold(ismember(unqfrq41, frequencies_local_optoOFF)) = frequencies_local_optoOFF;
            
            wait_times_all_optoOFF = [wait_times_all_optoOFF; wt_hold];
            frequencies_all_optoOFF = [frequencies_all_optoOFF; fl_hold];
            
        end
        
        % means and ses SMOOTHED
        mean_wait_times_optoON = nan(size(wait_times_all_optoON));
        std_wait_times_optoON = nan(size(mean_wait_times_optoON));
        se_wait_times_optoON = nan(size(mean_wait_times_optoON));
        mean_wait_times_optoOFF = nan(size(wait_times_all_optoOFF));
        std_wait_times_optoOFF = nan(size(mean_wait_times_optoOFF));
        se_wait_times_optoOFF = nan(size(mean_wait_times_optoOFF));
        
        for isesh = 1:size(wait_times_all_optoON,1)
            [mean_wait_times_optoON(isesh,:), std_wait_times_optoON(isesh,:), se_wait_times_optoON(isesh,:)] = nansmooth_adaptive_ampm(wait_times_all_optoON(isesh,:), 5, 3);
            [mean_wait_times_optoOFF(isesh,:), std_wait_times_optoOFF(isesh,:), se_wait_times_optoOFF(isesh,:)] = nansmooth_adaptive_ampm(wait_times_all_optoOFF(isesh,:), 5, 3);
        end

        % plot
        if plot_on >= 1
            hold on

            if plot_on==1
                if size(mean_wait_times_optoON,1)>1

                    %mean_wait_times_optoOFF= zscore_mtx(mean_wait_times_optoOFF')';
                    %mean_wait_times_optoON= zscore_mtx(mean_wait_times_optoON')';

                    plot(unqfrq41, mean(mean_wait_times_optoON) - std(mean_wait_times_optoON), 'linewidth', 1.0, 'color', colors(1,:))
                    plot(unqfrq41, mean(mean_wait_times_optoON) + std(mean_wait_times_optoON), 'linewidth', 1.0, 'color', colors(1,:))
                    plot(unqfrq41, mean(mean_wait_times_optoON) - std(mean_wait_times_optoON)./sqrt(size(mean_wait_times_optoON,1)), 'linewidth', 1.5, 'color', colors(1,:))
                    plot(unqfrq41, mean(mean_wait_times_optoON) + std(mean_wait_times_optoON)./sqrt(size(mean_wait_times_optoON,1)), 'linewidth', 1.5, 'color', colors(1,:))
                    plot(unqfrq41, mean(mean_wait_times_optoOFF) - std(mean_wait_times_optoOFF), 'linewidth', 1.0, 'color', colors(2,:))
                    plot(unqfrq41, mean(mean_wait_times_optoOFF) + std(mean_wait_times_optoOFF), 'linewidth', 1.0, 'color', colors(2,:))
                    plot(unqfrq41, mean(mean_wait_times_optoOFF) - std(mean_wait_times_optoOFF)./sqrt(size(mean_wait_times_optoOFF,1)), 'linewidth', 1.5, 'color', colors(2,:))
                    plot(unqfrq41, mean(mean_wait_times_optoOFF) + std(mean_wait_times_optoOFF)./sqrt(size(mean_wait_times_optoOFF,1)), 'linewidth', 1.5, 'color', colors(2,:))
                else
                    wait_times_plot_opto(trl_mtx,1);
                end

            elseif plot_on ==2
                if size(mean_wait_times_optoON,1)>1
                    errorbar(unqfrq41, mean(mean_wait_times_optoON), std(mean_wait_times_optoON)./sqrt(size(mean_wait_times_optoON,1)), 'color', colors(1,:))
                    errorbar(unqfrq41, mean(mean_wait_times_optoOFF), std(mean_wait_times_optoOFF)./sqrt(size(mean_wait_times_optoON,1)), 'color', colors(2,:))
                    %plot(unqfrq41, mean(mean_wait_times_optoON), 'linewidth', 1.5, 'color', colors(1,:))
                    %plot(unqfrq41, mean(mean_wait_times_optoON) - std(mean_wait_times_optoON)./sqrt(size(mean_wait_times_optoON,1)), 'linewidth', 1.0, 'color', colors(1,:))
                    %plot(unqfrq41, mean(mean_wait_times_optoON) + std(mean_wait_times_optoON)./sqrt(size(mean_wait_times_optoON,1)), 'linewidth', 1.0, 'color', colors(1,:))
                    
                else
                    wait_times_plot_opto(trl_mtx,1);
                end
            end
            %}


            %aesthetics
            set(gca,'TickLength',[0, 0]);
            xlabel('Tone Frequency (Hz)')
            ylabel('Wait Durations (s)')
            ylim_hold = ylim; ylim([0 ylim_hold(2)]);
            set(gca, 'XScale', 'log')
            title('Waits')
            ylim([0 40])
            xlim([4500 38500])
            xticks([5000 8500 14000 23000 35000])
        end
    end
    %}
end

if plot_on >= 1
    ylim([0 60])
    rich_bounds_2
    legend({'LaserON', 'LaserOFF'}, 'location', 'northeast')
end


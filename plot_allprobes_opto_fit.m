function [subj_curves, coefficients] = plot_allprobes_opto_smooth(file_paths, varargin)
% plot all preprobes

%get_file_paths_targeted('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\repeat_probes', {'preprobe', 'consistent'})


%% input prep
if nargin==2
    plot_on = varargin{1};
else
    plot_on = 1;
end

% number of fit coefficients
num_coefs = 4;

% smoothing instructions
smooth_min_bins = 5;
smooth_min_obs_pts = 3;

% opto tones
load('unqfrq41.mat', 'unqfrq41')
tones_optoON = unqfrq41(1:2:41);
tones_optoOFF = unqfrq41(2:2:40);

% points in fit curve
bins = 41;
fit_xvals = 1:length(unqfrq41)/bins:length(unqfrq41);

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


%% compute each probe
% iterate through unique probe numbers
for inum = max(pp_nums)
    
   % ppnum_fps = file_paths(pp_nums==inum); 
    %}
    ppnum_fps = file_paths;

    % gather wait times and fit
    %
    if ~isempty(ppnum_fps)
        
        %get wait times from every session
        wait_times_all_optoON = [];
        frequencies_all_optoON = [];
        wait_times_all_optoOFF = [];
        frequencies_all_optoOFF = [];
        coefficients = nan(size(ppnum_fps,1), num_coefs, 2); %3d is ON OFF
        subj_curves = nan(size(ppnum_fps,1), length(fit_xvals), 2);

        for ippnum = 1:size(ppnum_fps,1)
            
            % load
            load(ppnum_fps{ippnum}, 'trl_mtx')

            % compute means
            %

            % opto ON wait time
            [wait_times_local_optoON, frequencies_local_optoON] = wait_times_prep(trl_mtx(ismember(floor(trl_mtx(:,2)), tones_optoON),:), 1, 0);
            
            wt_hold = nan(1, length(unqfrq41));
            wt_hold(ismember(unqfrq41, frequencies_local_optoON)) = wait_times_local_optoON;
            fl_hold = nan(1, length(unqfrq41));
                    ns_hold = nansmooth_adaptive_ampm(wt_hold, smooth_min_bins, smooth_min_obs_pts);
                        while any(isnan(ns_hold))
                            ns_hold = nansmooth_adaptive_ampm(ns_hold, smooth_min_bins, smooth_min_obs_pts);
                        end
                    ns_hold(~isnan(wt_hold)) = wt_hold(~isnan(wt_hold));
                    wt_hold(1:2:41) = ns_hold(1:2:41);
                    fl_hold(1:2:41) = unqfrq41(1:2:41);

            wait_times_all_optoON = [wait_times_all_optoON; wt_hold];
            frequencies_all_optoON = [frequencies_all_optoON; fl_hold];

                % fit
                mean_wait_times = nanmean(wait_times_local_optoON);
                %try
                [~, coefficients(ippnum,:,1), modelFun_ON] = ampm_normal_logistic_fit_algo(find(~isnan(fl_hold)), wt_hold(~isnan(wt_hold)), [mean_wait_times 20 0 1]);                
                %catch
                %    figure; hold on
                %    plot(find(~isnan(wt_hold)), wt_hold(~isnan(wt_hold)), 'o')
                %    hold on; plot(xlim, [1 1].*mean_wait_times, 'k--')
                %    coefficients(ippnum,:,1) = [mean_wait_times zeros(1,num_coefs-1)];
                %end

                % santity check coefs
                if any(...
                        coefficients(ippnum,1,1) > 45 | coefficients(ippnum,1,1) < -5 | isnan(coefficients(ippnum,1,1)) | ...
                        coefficients(ippnum,2,1) > 41 | coefficients(ippnum,2,1) < -5 | isnan(coefficients(ippnum,2,1)) |...
                        coefficients(ippnum,3,1) > 25 | coefficients(ippnum,3,1) < -25 | isnan(coefficients(ippnum,3,1)) |...
                        coefficients(ippnum,4,1) > 25 | coefficients(ippnum,4,1) < -25 | isnan(coefficients(ippnum,4,1)) ...
                        )
                    disp(['Unreasonable optoON coefficients: ' num2str(coefficients(ippnum,:,1))])
                    coefficients(ippnum,:,1) = [mean_wait_times zeros(1,num_coefs-1)];
                end

                % load fit curve
                subj_curves(ippnum,:,1) = modelFun_ON(coefficients(ippnum,:,1), fit_xvals');


            % opto OFF
            [wait_times_local_optoOFF, frequencies_local_optoOFF] = wait_times_prep(trl_mtx(ismember(floor(trl_mtx(:,2)), tones_optoOFF),:), 1, 0);
            
            wt_hold = nan(1, length(unqfrq41));
            wt_hold(ismember(unqfrq41, frequencies_local_optoOFF)) = wait_times_local_optoOFF;
            fl_hold = nan(1, length(unqfrq41));
                    ns_hold = nansmooth_adaptive_ampm(wt_hold, smooth_min_bins, smooth_min_obs_pts);
                        while any(isnan(ns_hold))
                            ns_hold = nansmooth_adaptive_ampm(ns_hold, smooth_min_bins, smooth_min_obs_pts);
                        end
                    ns_hold(~isnan(wt_hold)) = wt_hold(~isnan(wt_hold));
                    wt_hold(2:2:40) = ns_hold(2:2:40);
                    fl_hold(2:2:40) = unqfrq41(2:2:40);

                % fit
                mean_wait_times = nanmean(wait_times_local_optoOFF);
                %try
                [~, coefficients(ippnum,:,2), modelFun_OFF] = ampm_normal_logistic_fit_algo(find(~isnan(fl_hold)), wt_hold(~isnan(wt_hold)), [mean_wait_times 20 0 1]);
                %catch
                %    coefficients(ippnum,:,2) = [mean_wait_times zeros(1,num_coefs-1)];
                %end
                    

                wait_times_all_optoOFF = [wait_times_all_optoOFF; wt_hold];
                frequencies_all_optoOFF = [frequencies_all_optoOFF; fl_hold];


                % santity check coefs
                if any(...
                        coefficients(ippnum,1,2) > 45 | coefficients(ippnum,1,2) < -5 | isnan(coefficients(ippnum,1,2)) | ...
                        coefficients(ippnum,2,2) > 41 | coefficients(ippnum,2,2) < -5 | isnan(coefficients(ippnum,2,2)) |...
                        coefficients(ippnum,3,2) > 25 | coefficients(ippnum,3,2) < -25 | isnan(coefficients(ippnum,3,2)) |...
                        coefficients(ippnum,4,2) > 25 | coefficients(ippnum,4,2) < -25 | isnan(coefficients(ippnum,4,2)) ...
                        )
                    disp(['Unreasonable optoOFF coefficients: ' num2str(coefficients(ippnum,:,2))])
                    coefficients(ippnum,:,2) = [mean_wait_times zeros(1,num_coefs-1)];
                end

                % load fit curve
                subj_curves(ippnum,:,2) = modelFun_OFF(coefficients(ippnum,:,2), fit_xvals');
                %plot(1:0.001:num_tones, subj_curves(iprobe,:,2), '-', 'color', colors(2,:), 'linewidth', 1)
                
        end

            % plot
            if plot_on==1
                figure; hold on
                if size(subj_curves,1)>1
                    plot(fit_xvals, mean(subj_curves(:,:,1)), 'linewidth', 2.0, 'color', colors(1,:))
                    plot(fit_xvals, mean(subj_curves(:,:,1)) - std(subj_curves(:,:,1))./sqrt(size(subj_curves(:,:,1),1)), 'linewidth', 1.0, 'color', colors(1,:))
                    plot(fit_xvals, mean(subj_curves(:,:,1)) + std(subj_curves(:,:,1))./sqrt(size(subj_curves(:,:,1),1)), 'linewidth', 1.0, 'color', colors(1,:))
                    plot(fit_xvals, mean(subj_curves(:,:,2)), 'linewidth', 2.0, 'color', colors(2,:))
                    plot(fit_xvals, mean(subj_curves(:,:,2)) - std(subj_curves(:,:,2))./sqrt(size(subj_curves(:,:,2),1)), 'linewidth', 1.0, 'color', colors(2,:))
                    plot(fit_xvals, mean(subj_curves(:,:,2)) + std(subj_curves(:,:,2))./sqrt(size(subj_curves(:,:,2),1)), 'linewidth', 1.0, 'color', colors(2,:))
                else
                    wait_times_plot_opto(trl_mtx,1);
                end

                %aesthetics
                set(gca,'TickLength',[0, 0]);
                xlabel('Tone Frequency (Hz)')
                ylabel('Wait Durations (s)')
                ylim_hold = ylim; ylim([0 ylim_hold(2)]);
                %set(gca, 'XScale', 'log')
                title('Waits')
                ylim([0 40])
                %xlim([4500 38500])
                %xticks([5000 8500 14000 23000 35000])
            end
    end
end
%}

if plot_on >= 1
    ylim([0 40])
    %rich_bounds_2
    legend({'LaserON', 'LaserOFF'}, 'location', 'northeast')
end


%% plot subject curves and dot plots
%
for i = 1:size(subj_curves,1)


    figure; hold on

    % plot laser ON
    plot(find(~isnan(wait_times_all_optoON(i,:))), wait_times_all_optoON(i, ~isnan(wait_times_all_optoON(i,:))), 'o', 'color', colors(1,:))
    plot(linspace(1, length(unqfrq41), length(subj_curves(i,:,1))), subj_curves(i,:,1), '-', 'color', colors(1,:))

    % plot laser OFF
    plot(find(~isnan(wait_times_all_optoOFF(i,:))), wait_times_all_optoOFF(i, ~isnan(wait_times_all_optoOFF(i,:))), 'o', 'color', colors(2,:))
    plot(linspace(1, length(unqfrq41), length(subj_curves(i,:,1))), subj_curves(i,:,2), '-', 'color', colors(2,:))

    title([num2str(i) '; ' num2str(coefficients(i,1,1)) '; ' num2str(coefficients(i,2,1)) '; ' num2str(coefficients(i,3,1)) '; ' num2str(coefficients(i,4,1))]);


end
%}

%% error bar of wait times
figure; hold on
wtON_mean = mean(wait_times_all_optoON(:,1:2:41),1);
wtON_se = std(wait_times_all_optoON(:,1:2:41), [], 1)./sqrt(size(wait_times_all_optoON,1));
errorbar(1:2:41, wtON_mean, wtON_se, 'color', colors(1,:), 'linewidth', 1)

wtOFF_mean = mean(wait_times_all_optoOFF(:,2:2:40), 1);
wtOFF_se = std(wait_times_all_optoOFF(:,2:2:40), [], 1)./sqrt(size(wait_times_all_optoOFF,1));
errorbar(2:2:40, wtOFF_mean, wtOFF_se, 'color', colors(2,:), 'linewidth', 1)

for i = 2:2:40
    [~,pval] = ttest(wait_times_all_optoOFF(:,i), mean([wait_times_all_optoON(:,i-1) wait_times_all_optoON(:,i+1)],2));
    sig_asterisks(pval, i, 35)
end

% find probe number
problem_num = str2num(file_paths{1}(strfind(file_paths{1}, 'opto_')+5:strfind(file_paths{1}, 'opto_')+6))-1;

ylim([0 40])
rich_bounds_prob('mevar', problem_num, 1);
set(gca,'TickLength',[0, 0]); box off;
xlim([-.9 42.9])


function [rvals_all] = tw_activity_corr_waits(trl_mtx, trl_idx, all_pkl_frames, traces, neurons, align_event, tw_bounds)
% runs tw_activity and produces one plot for each neuron comparing activity
% for each unique tone
%
% align events:
% NP onset = 1
% NP offset = 2
% Tone on = 3
% Head entry = 4
% Reward delivery = 5
% Reward receipt = 6
% Trial end = 7
%
% tw_bounds is a two item vector. the first item is the number of seconds
% before the align event the tw should start and the second item is the
% number of seconds after the align event the tw should end. eg [2 3] means
% the tw starts 2 seconds before and ends 3s after the align event.
% compute all activity

%activity matrix
trials = 1:size(trl_mtx,1);
act_mtx = tw_activity_II(trl_mtx, trl_idx, all_pkl_frames, traces, neurons, trials, align_event, tw_bounds);

% average across window
act_mtx_mean = nanmean(act_mtx,2);
act_mtx_mean = reshape(act_mtx_mean, size(act_mtx,1), size(act_mtx,3));

% probe trials only
trl_mtx_p = trl_mtx(trl_mtx(:,3)==0,:);
act_mtx_mean = act_mtx_mean(:,trl_mtx(:,3)==0);
act_mtx = act_mtx(:,:,trl_mtx(:,3)==0);
probe_waits = trl_mtx_p(:,12);

% compute r, p , and plot fit line for each cell
%
rpvals = nan(size(act_mtx,1), 2); 
for ic = 1:size(act_mtx_mean,1)
    
    % plot fit line
    figure; 
    [r, p] = fit_line(act_mtx_mean(ic,:)', probe_waits); 
    title(['neuron number ' num2str(ic) '; r=' num2str(r) ', p=' num2str(p)]);
    
    % aesthetics
    xlabel('mean activity over window')
    ylabel('wait time')
    
    % load values
    rpvals(ic,:) = [r p]; 
     
end
%}


% compute average rvalue at each time point over window
%

% preallocate rval matrix
rvals_all = nan(size(act_mtx,1), size(act_mtx,2));

% over all cells
for ic = 1:size(act_mtx,1)
    
    % over all time windows
    for iw = 1:size(act_mtx,2)

        % load values
        rvals_all(ic,iw) = corr(squeeze(act_mtx(ic,iw,:)), trl_mtx_p(:,12)); 
    
    end
end

rvals_all = abs(rvals_all);

% heatmap
%
figure; hold on;
imagesc(rvals_all)
plot(tw_bounds(1).*100.*[1 1], ylim, 'r-', 'linewidth', 2)
%}

% average
figure; 
hold on;
h = plot(nanfastsmooth(nanmean(rvals_all),11), 'linewidth', 2);
plot(nanfastsmooth(nanmean(rvals_all)+nanstd(rvals_all)./sqrt(size(rvals_all,1)),11), '-', 'color', get(h,'Color') , 'linewidth', 1)
plot(nanfastsmooth(nanmean(rvals_all)-nanstd(rvals_all)./sqrt(size(rvals_all,1)),11), '-', 'color', get(h,'Color') , 'linewidth', 1)
plot(tw_bounds(1).*100.*[1 1], ylim, 'r-', 'linewidth', 2)




    
function [all_deltas, all_waits, log_cell, norm_cell] = ALL_probe_over_probe_delta(training_group)
% plots probe to probe change in wait times

figure; 

% probe wait times
[~, ~, all_waits] = probe_wait_times(training_group, 1:8, 1:41);

% deltas
all_deltas = cell(1,length(all_waits)-1);
log_cell = cell(1,length(all_waits)-1);
norm_cell = cell(1,length(all_waits)-1);
for idelta = 1:7
    all_deltas{idelta} = all_waits{idelta+1} - all_waits{idelta};
end

% fit model to every delta
all_coefs = nan(size(all_deltas{idelta},1), 4, length(all_deltas));
for idelta = 1:7
    for isubj = 1:size(all_deltas{idelta},1)
        wait_deltas = all_deltas{idelta}(isubj,:);
        
        %{
        [~, all_coefs(isubj,:,idelta), ~] = ampm_normal_logistic_fit_algo(1:size(wait_deltas,2), wait_deltas);
        if any(squeeze(abs(all_coefs(isubj,1:4,idelta)))>30)
            all_coefs(isubj,1:4,idelta) = nan;
        end
        %}
        %
        try
            [~, all_coefs(isubj,:,idelta), modelFun] = ampm_normal_logistic_fit(1:size(wait_deltas,2), wait_deltas, [0 20 0 1]);
            %[~, all_coefs(isubj,1:3,idelta), modelFun] = ampm_logistic_fit(1:size(wait_deltas,2), wait_deltas, [0 20 0]);
            if any(squeeze(abs(all_coefs(isubj,1:4,idelta)))>30)
                all_coefs(isubj,1:4,idelta) = nan;
            end
        catch
            try
                %[~, all_coefs(isubj,1:3,idelta), modelFun] = ampm_logistic_fit(1:size(wait_deltas,2), wait_deltas, [mean(wait_deltas) 20.5 0]);
                [~, all_coefs(isubj,:,idelta), modelFun] = ampm_normal_logistic_fit(1:size(wait_deltas,2), wait_deltas, [mean(wait_deltas) 20.5 0 1]);
            catch
                all_coefs(isubj,1:4,idelta) = nan;
            end
        end
        %}
            
    end
end

% logistic
figure; 
subplot_idx = reshape(1:18, 3,6)';
for idelta = 1:6
    
    %heatmap
    subplot(6, 3, subplot_idx(idelta,1))
    imagesc(all_deltas{idelta})
    caxis([-15 15])
    colorbar
    if contains(training_group, 'hivar')
        rich_bounds_prob('hivar', idelta, 1);
    else
        rich_bounds_prob('mevar', idelta, 1);
    end

    % errorbar delta
    subplot(6, 3, subplot_idx(idelta,2))
    errorbar_mtx_lineonly(all_deltas{idelta})
    ylim([-11 11])
    hold on; plot(xlim, [0 0], 'k--')
    if contains(training_group, 'hivar')
        rich_bounds_prob('hivar', idelta, 1);
    else
        rich_bounds_prob('mevar', idelta, 1);
    end
    
    % errorbar_bar delta log
    if rem(idelta,2)>0 %if odd
        coefs_idelta = -squeeze(all_coefs(:,3,idelta));
    else
        coefs_idelta = squeeze(all_coefs(:,3,idelta));
    end
    log_cell{idelta} = coefs_idelta;

    
    % errorbar_bar delta norm
    norm_cell{idelta} = squeeze(all_coefs(:,4,idelta));
    
    % errorbar waits
    subplot(6, 3, subplot_idx(idelta,3))
    errorbar_mtx_lineonly(all_waits{idelta+1})
    ylim([0 25])
    hold on; plot(xlim, mean(mean(all_waits{1})).*[1 1], 'k--')
    if contains(training_group, 'hivar')
        rich_bounds_prob('hivar', idelta, 1);
    else
        rich_bounds_prob('mevar', idelta, 1);
    end

    
end

% preprobe waits
figure; hold on
errorbar_mtx_lineonly(all_waits{1})
hold on; plot(xlim, mean(mean(all_waits{1})).*[1 1], 'k--')
if contains(training_group, 'hivar')
    rich_bounds_prob('hivar', idelta, 1);
else
    rich_bounds_prob('mevar', idelta, 1);
end


% plot coef errorbars
figure; 
subplot(2,2, 1:2);
errorbar_plot(log_cell(1:6), 1); title log; hold on; plot(xlim, [0 0], 'k--')
subplot(2,2, 3:4);
errorbar_plot(norm_cell(1:6), 1); title norm; hold on; plot(xlim, [0 0], 'k--')


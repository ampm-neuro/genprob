function image_pca_plot_session(all_matrices, tses)
% plot all time bins colored by trial


%% compute activity matrix
%{
% load session first!!!
load('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_imaging_hpc\658648m2\LED_gen14_preprobe_06_01d.mat')
tses = [0.2 1.1 1.0 2.0 8.5 2.0];
[all_matrices, tones] = tw_activity_plot_tone_hm_warp_PCAversion(trl_mtx, trl_idx, frame_times, traces,...
    1:size(traces,1), 1:size(trl_mtx,1), 3, [5 25], tses);
%}

%% identify tones and trial stages

% combine dimensions
combo_tonebins = reshape(all_matrices, size(all_matrices,1)*size(all_matrices,2), size(all_matrices,3)); 



% build indices for plotting with color
tone_idx = repmat((1:size(all_matrices,1))', 1, size(all_matrices,2), size(all_matrices,3));
tone_idx = reshape(tone_idx, size(all_matrices,1)*size(all_matrices,2), size(all_matrices,3));

stage_idx = repmat((1:size(all_matrices,2)), size(all_matrices,1), 1, size(all_matrices,3));
stage_idx = reshape(stage_idx, size(all_matrices,1)*size(all_matrices,2), size(all_matrices,3));

% remove nan time bins (before and after trial)
tone_idx = tone_idx(~isnan(sum(combo_tonebins,2)), :);
stage_idx = stage_idx(~isnan(sum(combo_tonebins,2)), :); stage_idx = stage_idx-(min(stage_idx(~isnan(stage_idx)))-1);
combo_tonebins = combo_tonebins(~isnan(sum(combo_tonebins,2)), :);

% set stage index to 1:6 (stage numbers)
event_frame = cumsum(tses).*100;
for istage = 1:length(event_frame)
    stage_idx(stage_idx>=istage & stage_idx<=event_frame(istage)) = istage;
end

% pca
pcas = pca(zscore_mtx(combo_tonebins)');


figure; hold on; imagesc(pcas)
for irl = 1:length(event_frame)
    plot(xlim, event_frame(irl).*[1 1].*41, 'r-')
end


%% plot
colors = parula(41);

pcadims = [1 2 5];

figure; hold on; 
plot3(pcas(:,pcadims(1)), pcas(:,pcadims(2)), pcas(:,pcadims(3)), '.')
axes_hold = axis;
close

%
for istage = 1:length(unique(stage_idx(:,1)))
    figure; hold on; 
    
    for itone=1:length(unique(tone_idx(:,1)))
    
        
        % indices
        itone_idx = tone_idx(:,1)==itone;
        istage_idx = stage_idx(:,1)==istage;
        
        % plot transparent
        scatter1 = scatter3(pcas(itone_idx & ~istage_idx,pcadims(1)), pcas(itone_idx & ~istage_idx,pcadims(2)), pcas(itone_idx & ~istage_idx,pcadims(3)), 8, colors(itone,:), 'filled');  
        alpha(scatter1, 0.2)
        
        % plot stage
        scatter3(pcas(itone_idx & istage_idx,pcadims(1)), pcas(itone_idx & istage_idx,pcadims(2)), pcas(itone_idx & istage_idx,pcadims(3)), 8, colors(itone,:), 'filled'); 
        

    end
    
    axis(axes_hold)
    
    title(['Stage ' num2str(istage)])
    set(gca,'TickLength',[0, 0]); box off;
end
%}























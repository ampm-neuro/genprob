function image_pca_plot_session(all_matrices, tones, tses)
% plot all time bins colored by trial


%% compute activity matrix
%{
% load session first!!!
%load('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_imaging_hpc\658648m2\LED_gen14_preprobe_02_01d.mat')
%load('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_imaging_hpc\658648m2\LED_gen14_preprobe_04_01d.mat')
load('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_imaging_hpc\658648m2\LED_gen14_preprobe_06_01d.mat')
tses = [0.2 1.1 1.0 2.0 8.5 2.0];
[all_matrices, tones] = tw_activity_plot_tone_hm_warp_PCAversion(trl_mtx, trl_idx, frame_times, traces,...
    1:size(traces,1), 1:size(trl_mtx,1), 3, [5 25], tses);
%}



%% prepare input for dim red

% rotate to cell X time X trial
all_matrices = permute(all_matrices, [3,2,1]);

% remove empty time bins and empty trials
%tones = tones(sum(sum(~isnan(all_matrices)),2)>0);
all_matrices = all_matrices(:,sum(sum(~isnan(all_matrices)),3)>0, sum(sum(~isnan(all_matrices)),2)>0);

% remove time windows after end of random delay
event_frame = cumsum(tses).*100;
all_matrices = all_matrices(:, 1:event_frame(4) ,:);
size(all_matrices)

% zscore each cell
%
for ic = 1:size(all_matrices,1)
    hold_ic = all_matrices(ic,:,:);
    all_matrices(ic,:,:) = reshape(zscore_mtx(hold_ic(:)), [1, size(all_matrices,2), size(all_matrices,3)]);
end
%}



%% dim reduction
% be sure to add tensor toolbox to path
am_model = cp_als(tensor(all_matrices), 10); % fit CP model with 10 components
am_model = double(am_model);



%% stage index to 1:4 (Nose Poke, Movement, Head Entry, Tone On)
stage_idx = 1:size(am_model,2);
event_frame = cumsum(tses).*100;
for istage = 1:4
    stage_idx(stage_idx>=istage & stage_idx<=event_frame(istage)) = istage;
end




%% plot

dims = [1 2 3];
colors = parula(size(am_model,3));

%
for istage = 1:4
    figure; hold on; 
    
    for itrl=1:size(am_model,3)

        % indices
        istage_idx = stage_idx==istage;
        
        % plot transparent
        scatter1 = scatter3(am_model(dims(1), ~istage_idx, itrl), am_model(dims(2), ~istage_idx, itrl), am_model(dims(3), ~istage_idx, itrl), 8, colors(itrl,:), 'filled');  
        alpha(scatter1, 0.05)
        
        % plot stage
        scatter3(am_model(dims(1), istage_idx, itrl), am_model(dims(2), istage_idx, itrl), am_model(dims(3), istage_idx, itrl), 8, colors(itrl,:), 'filled'); 

    end
    
    grid on
    axis square
    title(['Stage ' num2str(istage)])
    set(gca,'TickLength',[0, 0]); box off;
end
%}























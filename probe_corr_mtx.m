function [ALL_probe_waits, corr_mtx, cm_17, offdiag_cell, comp2one_cell, comp2seven_cell] = probe_corr_mtx(training_group)
% correlation matrix of probes

num_freqs = 41;

green_color = [21 154 72]./255; blue_color =  [43 73 159]./255;
light_green_color = [162 201 129]./255; light_blue_color = [102 164 213]./255;
if contains(training_group, 'mevar')
   primary_color = green_color;
   secondary_color = light_green_color;
else
    primary_color = blue_color;
   secondary_color = light_blue_color;

end

ALL_probe_waits = cell(1,8);
ALL_subj_ids = cell(1,8);

for iprobe = 1:8
    figure
    [probe_waits, subjects_ids] = plot_probe_stages(iprobe, training_group); 
    title(['After problem ' num2str(iprobe-1)]); 
    ylim([0 30]); 

    ALL_probe_waits{iprobe} = reshape(probe_waits, num_freqs, length(probe_waits)/num_freqs)';
    ALL_subj_ids{iprobe} = unique(subjects_ids, 'rows' ,'stable');
end


% insert nans for missing subjects (uses first probe as master subj list)
for iprobe = 2:size(ALL_subj_ids,2)
    nan_hold = nan(size(ALL_subj_ids{1},1),num_freqs);
    nan_hold(ismember(ALL_subj_ids{1}, ALL_subj_ids{iprobe}, 'rows'), :) = ALL_probe_waits{iprobe};
    ALL_probe_waits{iprobe} = nan_hold;
    
end
for ifig = 1:8; close; end


% corr matrix for each subject
%
corr_mtx = nan(size(ALL_probe_waits,2), size(ALL_probe_waits,2), size(ALL_subj_ids{1}, 1));

for isubj = 1:size(ALL_subj_ids{1},1)
    for iprobe1 = 1:size(ALL_probe_waits,2)        
        for iprobe2 = 1:size(ALL_probe_waits,2)

            waits_1 = ALL_probe_waits{iprobe1}(isubj,:)';
            waits_2 = ALL_probe_waits{iprobe2}(isubj,:)';
            
            nnan_idx = ~isnan(waits_1) & ~isnan(waits_2);
            
            if isempty(waits_1(nnan_idx)) || isempty(waits_2(nnan_idx))
                continue
            end
            if sum(nnan_idx)>10
                corr_mtx(iprobe1, iprobe2, isubj) = corr(waits_1(nnan_idx), waits_2(nnan_idx));
            end
            
        end
    end
end
%}


figure; 
imagesc(nanmean(corr_mtx(1:7, 1:7, :),3))
axis square
set(gca,'TickLength',[0, 0]);
caxis([-.5 .5])
xticks(1:7); xticklabels(0:6)
yticks(1:7); yticklabels(0:6)
ylabel('Probe number')
xlabel('Probe number')
title([training_group '; average indiv subjs corrs'])
colorbar
%}

% probes 1-7
cm_17 = corr_mtx(1:7,1:7,:);

% errorbar plots
offdiag = nan(40,6); % off diagonal
comp2one = nan(40,6); % similarity of every probe to 1
comp2seven = nan(40,6); % similarity of every probe to 7
for isubj = 1:size(cm_17,3)
    
    % local subj
    cma = cm_17(:,:,isubj);
    
    % edges
    comp2one(isubj,:) = cma(1,2:7)';
    comp2seven(isubj,:) = cma(1:6,end);
    
    % diag
    cma(1,:) = []; cma(:,end) = [];
    offdiag(isubj,:) = cma(logical(eye(size(cma))));

end

% atanh
%{
offdiag = atanh(offdiag);
comp2one = atanh(comp2one);
comp2seven = atanh(comp2seven);
%}


%mat2cell
offdiag_cell = cell(1,6);
comp2one_cell = cell(1,6);
comp2seven_cell = cell(1,6);
for iprobe = 1:size(comp2one,2)
    offdiag_cell{iprobe} = offdiag(:,iprobe);
    comp2one_cell{iprobe} = comp2one(:,iprobe);
    comp2seven_cell{iprobe} = comp2seven(:,iprobe);
end

% errorbar plots
%{
% off diagonal
figure; hold on; 
errorbar_plot(offdiag_cell, 1, [], secondary_color, primary_color);
hold on; plot(xlim, [1 1].*0, 'k--')
xticklabels({'0&1', '1&2', '2&3', '3&4', '4&5', '5&6'}); 
xlabel('Probe comparison')
ylabel('Similarity (r)')
ylim([-1 1])
title('Compare adjacent probes')
axis square
%tic_vect = [ -.9 -.8 -.6 -.3 0 .3 .7 .9 .975 .995 ]; yticks(atanh(tic_vect)); yticklabels(tic_vect)

% similarity of every probe to 1 (test degree to which nearer-in-time 
% probes are similar to each other)
figure; hold on; 
errorbar_plot(comp2one_cell, 1, [], secondary_color, primary_color);
hold on; plot(xlim, [1 1].*0, 'k--')
xticks(1:6)
xlabel('Probe')
ylabel('Similarity (r)')
ylim([-1 1])
title('Compare to probe 1')
axis square
%tic_vect = [ -.9 -.8 -.6 -.3 0 .3 .7 .9 .975 .995 ]; yticks(atanh(tic_vect)); yticklabels(tic_vect)

% similarity of every probe to 7 (test degree to which responses converged
% over training)
figure; hold on; 
errorbar_plot(comp2seven_cell, 1, [], secondary_color, primary_color);
xlabel('Probe')
xticks(1:6)
xticklabels(0:5);
ylabel('Similarity (r)')
hold on; plot(xlim, [1 1].*0, 'k--')
ylim([-1 1])
title('Compare to probe 7')
axis square
%tic_vect = [ -.9 -.8 -.6 -.3 0 .3 .7 .9 .975 .995 ]; yticks(atanh(tic_vect)); yticklabels(tic_vect)
%}

%}




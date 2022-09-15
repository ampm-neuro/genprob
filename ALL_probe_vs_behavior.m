function ALL_probe_vs_behavior(training_group)
% Used to compare training problem behavior with probe behavior in
% every which way. produces many plots. additional ones found near the
% bottom


%% discrimination problem behavior
% problem behavior
[all_prefs, subj_idx, problem_idx, first_mtx, last_mtx, all_sesh_ct, ...
    learn_rate_mtx, delta_mtx, min_discrim_mtx, unq_subjects_beh, ...
    first_mtx_rich, first_mtx_poor, last_mtx_rich, last_mtx_poor] = ALL_cont_learn_delta_dot(0, training_group);
unq_subjects_beh = cell2mat(unq_subjects_beh);

[all_stage_learn_mtx] = ALL_stage_learn(['two_tone\' training_group], 1:6, 2);
first_day_d = squeeze(all_stage_learn_mtx(:,1,:));
last_day_d = squeeze(all_stage_learn_mtx(:,2,:));

%% problem tones
load('unqfrq41', 'unqfrq41')
all_prob_tones = rich_bounds_prob('mevar', 0);
for iprob = 1:size(all_prob_tones,1)
   tones_hold = find(ismember(unqfrq41,all_prob_tones(iprob,:)));
   all_prob_tones(iprob,1) = tones_hold(ismember(tones_hold, 16:26));
   all_prob_tones(iprob,2) = tones_hold(ismember(tones_hold, setdiff(1:41, 16:26))); 
end



%% probe behavior

%
% probe
last_probe = 8;
rich_tones = 16:26; poor_tones = setdiff(1:41, rich_tones); 
figure; [rich_waits, unq_subjects_probe_rich] = probe_wait_times_nansmooth(training_group, 1:last_probe, rich_tones, 0);
[poor_waits] = probe_wait_times(training_group, 1:last_probe, poor_tones, 0);
[~, unq_subjects_probe_all, all_waits] = probe_wait_times(training_group, 1:last_probe, 1:41); % probe delta prep CHOOSE BEHAVIOR FOLDER
%close

for i = 1:length(unq_subjects_probe_all)
    
    % rich tone preference during probe tests
    diff_waits{i} = (rich_waits{i}-poor_waits{i})./(rich_waits{i}+poor_waits{i}); 
    dwaits_hold = nan(size(unq_subjects_beh,1),1);   
    dwaits_hold(ismember(unq_subjects_beh, unq_subjects_probe_rich{i}, 'rows')) = diff_waits{i}(ismember(unq_subjects_probe_rich{i}, unq_subjects_beh, 'rows'), :);
    diff_waits{i} = dwaits_hold;

    % probe to probe change
    if ismember(i, 1:last_probe-1)
        
        if length(unq_subjects_probe_all) < i+1
            continue
        end
        
        %[unq_subj_probe_delta{i}, pre_idx, post_idx] = intersect(unq_subjects_probe_all{i},unq_subjects_probe_all{i+1},'rows', 'stable');
        
        % raw data (should contains NaNs - no extrapolation)
        probes_pre = all_waits{i};%(pre_idx,:);
        probes_post = all_waits{i+1};%(post_idx,:);

        
        % plot both probe wait times (raw)
        figure; hold on
        errorbar(nanmean(probes_pre,1), nanstd(probes_pre,[],1)./sqrt(sum(~isnan(probes_pre),1)));
        errorbar(nanmean(probes_post,1), nanstd(probes_post,[],1)./sqrt(sum(~isnan(probes_post),1)));
        ylim([0 40])
        title(['Probes ' num2str(i) ' and ' num2str(i+1)])

        % smooth probe wait times
        %
        for ismooth = 1:size(probes_pre,1)
            probes_pre(ismooth,:) = nansmooth_ampm(probes_pre(ismooth,:), 11); %removes nans
            probes_post(ismooth,:) = nansmooth_ampm(probes_post(ismooth,:), 11);
        end
        %}

        % similarity measures(correlation and difference)
        %
        probe_delta_r{i} = nan(size(probes_pre,1),1); % correlation
        
        for icorr = 1:size(probes_pre,1)
            
            % check missing data
            nnan_idx = ~isnan(probes_pre(icorr,:)) & ~isnan(probes_post(icorr,:));
            if sum(nnan_idx)==0
                probe_delta_r{i}(icorr) = nan;
                continue
            end
            
            probe_delta_r{i}(icorr) = corr(probes_pre(icorr, nnan_idx)', probes_post(icorr,nnan_idx)');
        end
        
        %
        pd_mtx = probes_post-probes_pre; %normalized difference
        probe_delta_norm{i} = nansum(abs(pd_mtx),2) ./ nansum(abs(probes_pre)+abs(probes_post),2);


        % plot probe change surrounding each problem
        %
        pd_mtx = pd_mtx./(probes_post+probes_pre);

        % pre probes only
        if i < 7
        
            % change outside the poor and rich tones
            if all_prob_tones(i,1)>all_prob_tones(i,2)
                bpc_idx = 1:all_prob_tones(i,2);
                brc_idx = all_prob_tones(i,1):size(pd_mtx,2);
            else
                bpc_idx = all_prob_tones(i,2):size(pd_mtx,2);
                brc_idx = 1:all_prob_tones(i,1);
            end
                beyond_poor_change = mean(pd_mtx(:,bpc_idx), 2);
                beyond_rich_change = mean(pd_mtx(:,brc_idx), 2);
            
        end

            
        % plot probe over probe change in smoothed normalized wait times
        figure; hold on;
        plot(nanmean((pd_mtx)), 'k-', 'linewidth', 2)
        plot(nanmean((pd_mtx)) + nanstd((pd_mtx))./sqrt(sum(~isnan(pd_mtx))), 'k-', 'linewidth', 1);
        plot(nanmean((pd_mtx)) - nanstd((pd_mtx))./sqrt(sum(~isnan(pd_mtx))), 'k-', 'linewidth', 1);
        ylim([-.25 .25])
        xlim([0 42])
        title(['Change from probes ' num2str(i) ' to ' num2str(i+1)])
        hold on; plot(xlim, [0 0], 'k--')
        % pre probes only
        if i < 7
            hold on; plot(all_prob_tones(i,1).*[1 1], ylim, 'r-')
            hold on; plot(all_prob_tones(i,2).*[1 1], ylim, 'k-')
        end
        set(gca,'TickLength',[0, 0]); box off;

        % only subjects that have relevant problem data
        %{
        probe_delta_hold_r = nan(size(unq_subjects_beh,1),1);
        probe_delta_hold_r(ismember(unq_subjects_beh, unq_subj_probe_delta{i}, 'rows')) = probe_delta_r{i}(ismember(unq_subj_probe_delta{i}, unq_subjects_beh, 'rows'), :);
        probe_delta_r{i} = probe_delta_hold_r;
        
        probe_delta_hold_norm = nan(size(unq_subjects_beh,1),1);
        probe_delta_hold_norm(ismember(unq_subjects_beh, unq_subj_probe_delta{i}, 'rows')) = probe_delta_norm{i}(ismember(unq_subj_probe_delta{i}, unq_subjects_beh, 'rows'), :);
        probe_delta_norm{i} = probe_delta_hold_norm;

        probe_delta_br_hold = nan(size(unq_subjects_beh,1),1);
        probe_delta_bp_hold = nan(size(unq_subjects_beh,1),1);
        probe_delta_br_hold(ismember(unq_subjects_beh, unq_subj_probe_delta{i}, 'rows')) = beyond_rich_change(ismember(unq_subj_probe_delta{i}, unq_subjects_beh, 'rows'));
        probe_delta_bp_hold(ismember(unq_subjects_beh, unq_subj_probe_delta{i}, 'rows')) = beyond_poor_change(ismember(unq_subj_probe_delta{i}, unq_subjects_beh, 'rows'));
        probe_delta_br{i} = probe_delta_br_hold;
        probe_delta_bp{i} = probe_delta_bp_hold;
        %}
        probe_delta_br{i} = beyond_rich_change;
        probe_delta_bp{i} = beyond_poor_change;
    end

end

diff_waits = cell2mat(diff_waits)
diff_waits_past = diff_waits(:,1:length(unq_subjects_probe_all)-1);
diff_waits_future = diff_waits(:,2:length(unq_subjects_probe_all));
probe_delta_r = cell2mat(probe_delta_r);
probe_delta_norm = cell2mat(probe_delta_norm);

% plot mean probe over probe change in similarlity (correlation)
% pd_full = probe_delta_r(~isnan(sum(probe_delta_r,2)),:); % complete only
pd_full = probe_delta_r;
pd_full = atanh(pd_full); %fischer z transform for normality
pd_full_vect = [];
for i = 1:length(unq_subjects_probe_all)-1
    pd_full_vect = [pd_full_vect {pd_full(:,i)}];
end
figure; errorbar_plot(pd_full_vect, 1)
title('mean probe over probe change')
xlabel('Training problem number')
ylabel('Probe over probe similarity (r)')
hold on; plot(xlim, [0 0], 'k--')
ylim([-1 1])
maxr = .975; tic_vect = [-maxr -.9 -.8 -.6 -.3 0 .3 .6 .8 .9 maxr];
ylim(atanh([-maxr maxr])); yticks(atanh(tic_vect)); yticklabels(tic_vect)
%[~, b, c, d] = ttest(mean(pd_full(:,[1 2 3]),2), mean(pd_full(:,[4 5 6]),2))

figure; errorbar_barplot([pd_full_vect(1) pd_full_vect(end)], 1);
%figure; errorbar_barplot([{atanh(pd_full_vect{1})} {atanh(pd_full_vect{end})}]);
title('Similarity: first two and last two probes')
ylabel('Probe over probe similarity (r)')
ylim([-1 1])
[~, r_FirstLast] = ttest(pd_full_vect{1}, pd_full_vect{end})
ylim(atanh([-maxr maxr])); yticks(atanh(tic_vect)); yticklabels(tic_vect)
%[~, r_FirstLast_z] = ttest(atanh(pd_full_vect{1}), atanh(pd_full_vect{end}))
% outside rich change greater than outside poor change
try
    probe_delta_br = cell2mat(probe_delta_br);
    probe_delta_bp = cell2mat(probe_delta_bp);
catch
    probe_delta_br_hold = [];
    probe_delta_bp_hold = [];
    for i1 = 1:length(probe_delta_br)
        probe_delta_br_hold = [probe_delta_br_hold; probe_delta_br{i1}];
        probe_delta_bp_hold = [probe_delta_bp_hold; probe_delta_bp{i1}];
    end
    probe_delta_br = probe_delta_br_hold;
    probe_delta_bp = probe_delta_bp_hold;
end
    % line errorbar plot
    title('rich side and poor side changes over')
    all_rich_changes = [];
    all_poor_changes = [];
    for iproblem = 1:length(unq_subjects_probe_all)-1
        all_rich_changes = [all_rich_changes {probe_delta_br(:,iproblem)}];
        all_poor_changes = [all_poor_changes {probe_delta_bp(:,iproblem)}];
    end
    figure; hold on;
    errorbar_plot(all_rich_changes(1:6), 1, [], [ 0                     0.447                     0.741])
    errorbar_plot(all_poor_changes(1:6), 1, [], [0.85                     0.325                     0.098])
    hold on; plot(xlim, [1 1].*0, 'k--')
        
colors = bone(10); colors = colors(2:last_probe,:);
figure; hold on;
for iproblem = 1:length(unq_subjects_probe_all)-1
    errorbar_plot([{probe_delta_bp(:,iproblem)}, {probe_delta_br(:,iproblem)}], 1, [1 2], colors(iproblem,:));
end
title('colors for errorbar plot')
figure
errorbar_barplot([{probe_delta_bp(:)}, {probe_delta_br(:)}], 1);
plot(xlim, [0 0], 'k--')
xticks([1 2])
xticklabels({'poorside', 'richside'})
title('rich change greater than poor change')
%[~,br1,br2,br3] = ttest(probe_delta_br(:))
%[~,bp1,bp2,bp3] = ttest(probe_delta_bp(:))



%% does change correlate with distance (in SD) between problem and probe responses

% problem responses output above
%[first_mtx_rich, first_mtx_poor, last_mtx_rich, last_mtx_poor]

% preallocate distances
dist_rich_preprobe_to_earlydiscr = nan(size(first_mtx)); % 17 x 6
dist_poor_preprobe_to_earlydiscr = nan(size(dist_rich_preprobe_to_earlydiscr));
dist_rich_preprobe_to_latediscr  = nan(size(dist_rich_preprobe_to_earlydiscr));
dist_poor_preprobe_to_latediscr = nan(size(dist_rich_preprobe_to_earlydiscr));
dist_rich_postprobe_to_earlydiscr = nan(size(dist_rich_preprobe_to_earlydiscr));
dist_poor_postprobe_to_earlydiscr = nan(size(dist_rich_preprobe_to_earlydiscr));
dist_rich_postprobe_to_latediscr  = nan(size(dist_rich_preprobe_to_earlydiscr));
dist_poor_postprobe_to_latediscr = nan(size(dist_rich_preprobe_to_earlydiscr));

% stdev between problem wait times and smoothed probe response mean at relevant tone
for iprob = 1:length(unq_subjects_probe_all)-1
    
    % tone indices
    rich_tone_idx = all_prob_tones(1);
    poor_tone_idx = all_prob_tones(2);
    
    % prep probe response data
    %
        [unq_subj_probe_delta{i}, pre_idx, post_idx]...
            = intersect(unq_subjects_probe_all{iprob},unq_subjects_probe_all{iprob+1},'rows', 'stable');

        % raw wait data, just relevant subjs (should contains NaNs - no extrapolation)
        probes_pre = all_waits{iprob}(pre_idx,:);
            probes_pre_std = nan(size(probes_pre));
        probes_post = all_waits{iprob+1}(post_idx,:);
            probes_post_std = nan(size(probes_post));

        % smooth each subj's waits
        for isubj = 1:size(probes_pre,1)
            [probes_pre(isubj,:), probes_pre_std(isubj,:)] = nansmooth_ampm(probes_pre(isubj,:), 11); %removes nans
            [probes_post(isubj,:), probes_post_std(isubj,:)] = nansmooth_ampm(probes_post(isubj,:), 11);
        end
    
    % problem subject index
    [~,problem_subj_idx, probe_subj_idx] = intersect(unq_subjects_beh, unq_subj_probe_delta{i},'rows', 'stable');

    % pre probe responses at that tone
    probe_pre_rich_mean = probes_pre(probe_subj_idx, rich_tone_idx);
    probe_pre_rich_std = probes_pre_std(probe_subj_idx, rich_tone_idx);
    probe_pre_poor_mean = probes_pre(probe_subj_idx, poor_tone_idx);
    probe_pre_poor_std = probes_pre_std(probe_subj_idx, poor_tone_idx);

    % post probe responses at that tone
    probe_post_rich_mean = probes_post(probe_subj_idx, rich_tone_idx);
    probe_post_rich_std = probes_post_std(probe_subj_idx, rich_tone_idx);
    probe_post_poor_mean = probes_post(probe_subj_idx, poor_tone_idx);
    probe_post_poor_std = probes_post_std(probe_subj_idx, poor_tone_idx);

    if iprob > 6
        continue
    end
    
    % compute and load distances
    dist_rich_preprobe_to_earlydiscr(problem_subj_idx,iprob) = ...
        (first_mtx_rich(problem_subj_idx,iprob)-probe_pre_rich_mean)./probe_pre_rich_std;
    dist_poor_preprobe_to_earlydiscr(problem_subj_idx,iprob) = ...
        (first_mtx_poor(problem_subj_idx,iprob)-probe_pre_poor_mean)./probe_pre_poor_std;
    dist_rich_preprobe_to_latediscr(problem_subj_idx,iprob) = ...
        (last_mtx_rich(problem_subj_idx,iprob)-probe_pre_rich_mean)./probe_pre_rich_std;
    dist_poor_preprobe_to_latediscr(problem_subj_idx,iprob) = ...
        (last_mtx_poor(problem_subj_idx,iprob)-probe_pre_poor_mean)./probe_pre_poor_std;
    
    dist_rich_postprobe_to_earlydiscr(problem_subj_idx,iprob) = ...
        (first_mtx_rich(problem_subj_idx,iprob)-probe_post_rich_mean)./probe_post_rich_std;
    dist_poor_postprobe_to_earlydiscr(problem_subj_idx,iprob) = ...
        (first_mtx_poor(problem_subj_idx,iprob)-probe_post_poor_mean)./probe_post_poor_std;
    dist_rich_postprobe_to_latediscr(problem_subj_idx,iprob) = ...
        (last_mtx_rich(problem_subj_idx,iprob)-probe_post_rich_mean)./probe_post_rich_std;
    dist_poor_postprobe_to_latediscr(problem_subj_idx,iprob) = ...
        (last_mtx_poor(problem_subj_idx,iprob)-probe_post_poor_mean)./probe_post_poor_std;


end



%% plots and stats

% colors
colors = bone(10); colors = colors(2:7,:);

% mixed model prep
subj_num_mtx = repmat((1:size(diff_waits_past,1))', 1, size(diff_waits_past,2)); % subject matrix\ 
    subj_num_mtx = subj_num_mtx(:,1:end-1);
stage_num_mtx = repmat(1:size(diff_waits_past,2), size(diff_waits_past,1), 1); % stage matrix
    stage_num_mtx = stage_num_mtx(:,1:end-1);
model_str = 'dv~ProbeBeh+(1|subject)+(ProbeBeh-1|subject)+(1|problem)+(ProbeBeh-1|problem)';
%model_str = 'dv~ProbeBeh+(1|subject)+(1|problem)';


%% probe changes as a function of problem-probe dissagreements

% distance between late rich problem mean and pre-probe VS postprobe rich
%

    % mean centered correlation dot plot
    figure; hold on; 
    
    var1 = dist_rich_preprobe_to_latediscr;
    var2 = probe_delta_br;
    for iprob = 1:length(unq_subjects_probe_all)-2
        % mean center
        %
        var1(:, iprob) = var1(:, iprob) - nanmean(var1(:, iprob));
        var2(:, iprob) = var2(:, iprob) - nanmean(var2(:, iprob));
        %}

        %[p,r] = corr(var1(:, iprob), var2(:, iprob))
        
        % old
        
        size(var1)
        size(var2)
        
        plot(var1(:, iprob), var2(:, iprob), 'o', 'color', colors(iprob,:))
        %plot(nanmean(var1(:, iprob)), nanmean(var2(:, iprob)), '.', 'color', colors(iprob,:), 'markersize', 50)
    end
    fitline_hold_mc = var1(:);
    probe_hold_mc = var2(:,1:end-1);
    fitline_hold_mc = fitline_hold_mc(1:length(probe_hold_mc(:)));
    fitline_hold = dist_rich_preprobe_to_latediscr(:);
    probe_hold = probe_delta_br(:,1:end-1);
    fitline_hold = fitline_hold(1:length(probe_hold(:)));
    [r, p] = fit_line(fitline_hold_mc, probe_hold_mc(:));
    set(gca,'TickLength',[0, 0]); box off;
    title(['Dissagreement and change (rich); r=' num2str(r) ', p=' num2str(p)])
    xlabel(['rich-side dissagreement (preprobe and late discrim)'])
    ylabel(['rich-side probe change (from pre to post)'])
    hold on; plot(xlim, [0 0], 'k--')
    hold on; plot([0 0], ylim, 'k--')
    %}
    
    
    % change over problem set dot plot
    figure; hold on;
    for isubj = 1:size(dist_rich_preprobe_to_latediscr,1)
        for iprob = 1:length(unq_subjects_probe_all)-2
            if iprob < length(unq_subjects_probe_all)-2
                plot(dist_rich_preprobe_to_latediscr(isubj, [iprob iprob+1]), probe_delta_br(isubj, [iprob iprob+1]), '-', 'color', colors(iprob,:), 'linewidth', 0.5)
            end
            plot(dist_rich_preprobe_to_latediscr(isubj, iprob), probe_delta_br(isubj, iprob), 'o', 'color', colors(iprob,:))
        end
    end
    r_and_p = nan(length(unq_subjects_probe_all)-2, 2);
    for iprob = 1:length(unq_subjects_probe_all)-2
        if iprob < length(unq_subjects_probe_all)-2
            plot(nanmean(dist_rich_preprobe_to_latediscr(:, [iprob iprob+1])), nanmean(probe_delta_br(:, [iprob iprob+1])), '-', 'color', colors(iprob,:), 'linewidth', 6)
        end
        plot(nanmean(dist_rich_preprobe_to_latediscr(:, iprob)), nanmean(probe_delta_br(:, iprob)), '.', 'color', colors(iprob,:), 'markersize', 60)
        [r_and_p(iprob,1), r_and_p(iprob,2)] = fit_line(dist_rich_preprobe_to_latediscr(:, iprob), probe_delta_br(:, iprob));
    end
    set(gca,'TickLength',[0, 0]); box off;
    title(['Dissagreement and change (rich); r=' num2str(r) ', p=' num2str(p)])
    xlabel(['rich-side dissagreement (preprobe and late discrim)'])
    ylabel(['rich-side probe change (from pre to post)'])
    hold on; plot(xlim, [0 0], 'k--')
    hold on; plot([0 0], ylim, 'k--')

            % mixed model
            tbl = table(probe_hold(:), fitline_hold, subj_num_mtx(:), stage_num_mtx(:), 'VariableNames',{'dv','ProbeBeh','subject','problem'});
            tbl.subject = categorical(tbl.subject);
            tbl.problem = categorical(tbl.problem);
            lme_dissagreement_rich = fitlme(tbl, model_str)
            
    % r val bar     
     figure; hold on
     bar(r_and_p(:,1)); 
     sig_asterisks(r_and_p(:,2), 1:length(r_and_p(:,2)), r_and_p(:,1).*1.2)
     ylim([-1.09 1.09])
     title(['dissagreement rich'])
        
        
% distance between late rich problem mean and pre-probe VS postprobe rich
%

        % mean centered correlation dot plot
        figure; hold on; 
        var1 = dist_poor_preprobe_to_latediscr;
        var2 = probe_delta_bp;
        for iprob = 1:length(unq_subjects_probe_all)-2
            % mean center
            %
            var1(:, iprob) = var1(:, iprob) - nanmean(var1(:, iprob));
            var2(:, iprob) = var2(:, iprob) - nanmean(var2(:, iprob));
            %}

            %[p,r] = corr(var1(:, iprob), var2(:, iprob))

            % old
            plot(var1(:, iprob), var2(:, iprob), 'o', 'color', colors(iprob,:))
            %plot(nanmean(var1(:, iprob)), nanmean(var2(:, iprob)), '.', 'color', colors(iprob,:), 'markersize', 50)
        end
        fitline_hold_mc = var1(:);
        probe_hold_mc = var2(:,1:end-1);
        fitline_hold_mc = fitline_hold_mc(1:length(probe_hold_mc(:)));
        fitline_hold = dist_poor_preprobe_to_latediscr(:);
        probe_hold = probe_delta_bp(:,1:end-1);
        fitline_hold = fitline_hold(1:length(probe_hold(:)));
        [r, p] = fit_line(fitline_hold_mc, probe_hold_mc(:));
        set(gca,'TickLength',[0, 0]); box off;
        title(['Dissagreement and change (poor); r=' num2str(r) ', p=' num2str(p)])
        xlabel(['poor-side dissagreement (preprobe and late discrim)'])
        ylabel(['poor-side probe change (from pre to post)'])
        hold on; plot(xlim, [0 0], 'k--')
        hold on; plot([0 0], ylim, 'k--')


    % change over problem set dot plot
    figure; hold on;
    for isubj = 1:size(dist_poor_preprobe_to_latediscr,1)
        for iprob = 1:length(unq_subjects_probe_all)-2
            if iprob < length(unq_subjects_probe_all)-2
                plot(dist_poor_preprobe_to_latediscr(isubj, [iprob iprob+1]), probe_delta_bp(isubj, [iprob iprob+1]), '-', 'color', colors(iprob,:), 'linewidth', 0.5)
            end
            plot(dist_poor_preprobe_to_latediscr(isubj, iprob), probe_delta_bp(isubj, iprob), 'o', 'color', colors(iprob,:))
        end
    end
    r_and_p = nan(length(unq_subjects_probe_all)-2, 2);
    for iprob = 1:length(unq_subjects_probe_all)-2
        if iprob < length(unq_subjects_probe_all)-2
            plot(nanmean(dist_poor_preprobe_to_latediscr(:, [iprob iprob+1])), nanmean(probe_delta_bp(:, [iprob iprob+1])), '-', 'color', colors(iprob,:), 'linewidth', 6)
        end
        plot(nanmean(dist_poor_preprobe_to_latediscr(:, iprob)), nanmean(probe_delta_bp(:, iprob)), '.', 'color', colors(iprob,:), 'markersize', 60)
        [r_and_p(iprob,1), r_and_p(iprob,2)] = fit_line(dist_poor_preprobe_to_latediscr(:, iprob), probe_delta_bp(:, iprob));
    end
    set(gca,'TickLength',[0, 0]); box off;
    title(['Dissagreement and change (poor); r=' num2str(r) ', p=' num2str(p)])
    xlabel(['poor-side dissagreement (preprobe and late discrim)'])
    ylabel(['poor-side probe change (from pre to post)'])
    hold on; plot(xlim, [0 0], 'k--')
    hold on; plot([0 0], ylim, 'k--')

            % mixed model
            tbl = table(probe_hold(:), fitline_hold, subj_num_mtx(:), stage_num_mtx(:), 'VariableNames',{'dv','ProbeBeh','subject','problem'});
            tbl.subject = categorical(tbl.subject);
            tbl.problem = categorical(tbl.problem);
            lme_dissagreement_poor = fitlme(tbl, model_str)
            
     % r val bar     
     figure; hold on
     bar(r_and_p(:,1)); 
     sig_asterisks(r_and_p(:,2), 1:length(r_and_p(:,2)), r_and_p(:,1).*1.2)
     ylim([-1.09 1.09])
     title(['dissagreement poor'])
        

%% PAST PROBE plots


% number of trials to crit
%
        % mean centered correlation dot plot
        figure; hold on; 
        var1 = diff_waits_past;
        var2 = all_sesh_ct;
        for iprob = 1:length(unq_subjects_probe_all)-2
            % mean center
            %
            var1(:, iprob) = var1(:, iprob) - nanmean(var1(:, iprob));
            var2(:, iprob) = var2(:, iprob) - nanmean(var2(:, iprob));
            %}

            %[p,r] = corr(var1(:, iprob), var2(:, iprob))

            % old
            plot(var1(:, iprob), var2(:, iprob), 'o', 'color', colors(iprob,:))
            %plot(nanmean(var1(:, iprob)), nanmean(var2(:, iprob)), '.', 'color', colors(iprob,:), 'markersize', 50)
        end
        fitline_hold_mc = var1(:);
        probe_hold_mc = var2(:,1:end-1);
        fitline_hold_mc = fitline_hold_mc(1:length(probe_hold_mc(:)));
        fitline_hold = all_sesh_ct(:);
        probe_hold = diff_waits_past(:,1:end-1);
        fitline_hold = fitline_hold(1:length(probe_hold(:)));
        [r, p] = fit_line(fitline_hold_mc, probe_hold_mc(:));
        set(gca,'TickLength',[0, 0]); box off;
        title(['num trials to crit; r=' num2str(r) ', p=' num2str(p)])
        xlabel(['prior probe discrim'])
        ylabel(['num trials to crit'])
        hold on; plot(xlim, [0 0], 'k--')
        hold on; plot([0 0], ylim, 'k--')

        
        % change over problem set dot plot
        figure; hold on;
        for isubj = 1:size(diff_waits_past,1)
            for iprob = 1:length(unq_subjects_probe_all)-2
                if iprob < length(unq_subjects_probe_all)-2
                    plot(diff_waits_past(isubj, [iprob iprob+1]), all_sesh_ct(isubj, [iprob iprob+1]), '-', 'color', colors(iprob,:), 'linewidth', 0.5)
                end
                plot(diff_waits_past(isubj, iprob), all_sesh_ct(isubj, iprob), 'o', 'color', colors(iprob,:))
            end
        end
        r_and_p = nan(length(unq_subjects_probe_all)-2, 2);%all r vals and p vals
        for iprob = 1:length(unq_subjects_probe_all)-2
            if iprob < length(unq_subjects_probe_all)-2
                plot(nanmean(diff_waits_past(:, [iprob iprob+1])), nanmean(all_sesh_ct(:, [iprob iprob+1])), '-', 'color', colors(iprob,:), 'linewidth', 6)
            end
            plot(nanmean(diff_waits_past(:, iprob)), nanmean(all_sesh_ct(:, iprob)), '.', 'color', colors(iprob,:), 'markersize', 60)
            [r_and_p(iprob,1), r_and_p(iprob,2)] = fit_line(diff_waits_past(:, iprob), all_sesh_ct(:, iprob));
        end
        set(gca,'TickLength',[0, 0]); box off;
        title(['num trials to crit'])
        xlabel(['prior probe discrim'])
        ylabel(['num trials to crit'])
        hold on; plot(xlim, [0 0], 'k--')
        hold on; plot([0 0], ylim, 'k--')

                % mixed model
                tbl = table(probe_hold(:), fitline_hold, subj_num_mtx(:), stage_num_mtx(:), 'VariableNames',{'dv','ProbeBeh','subject','problem'});
                tbl.subject = categorical(tbl.subject);
                tbl.problem = categorical(tbl.problem);
                lme_trials_to_crit = fitlme(tbl, model_str)
                
        % r val bar     
        figure; hold on
        bar(r_and_p(:,1)); 
        sig_asterisks(r_and_p(:,2), 1:length(r_and_p(:,2)), r_and_p(:,1).*1.2)
        ylim([-1.09 1.09])
        title(['num trials to crit'])
      

% early discrim (proportion)
%{
        % mean centered correlation dot plot
        figure; hold on; 
        var1 = diff_waits_past;
        var2 = first_mtx;
        for iprob = 1:length(unq_subjects_probe_all)-2
            % mean center
            %
            var1(:, iprob) = var1(:, iprob) - nanmean(var1(:, iprob));
            var2(:, iprob) = var2(:, iprob) - nanmean(var2(:, iprob));
            %

            %[p,r] = corr(var1(:, iprob), var2(:, iprob))

            % old
            plot(var1(:, iprob), var2(:, iprob), 'o', 'color', colors(iprob,:))
            %plot(nanmean(var1(:, iprob)), nanmean(var2(:, iprob)), '.', 'color', colors(iprob,:), 'markersize', 50)
        end
        fitline_hold_mc = var1(:);
        probe_hold_mc = var2(:,1:end-1);
        fitline_hold_mc = fitline_hold_mc(1:length(probe_hold_mc(:)));
        fitline_hold = first_mtx(:);
        probe_hold = diff_waits_past(:,1:end-1);
        fitline_hold = fitline_hold(1:length(probe_hold(:)));
        [r, p] = fit_line(fitline_hold_mc, probe_hold_mc(:));
        set(gca,'TickLength',[0, 0]); box off;
        title(['early discrim (proportion); r=' num2str(r) ', p=' num2str(p)])
        xlabel(['prior probe discrim'])
        ylabel(['early discrim (proportion)'])
        hold on; plot(xlim, [0 0], 'k--')
        hold on; plot([0 0], ylim, 'k--')
        
        
        % change over problem set dot plot
        figure; hold on;
        for isubj = 1:size(diff_waits_past,1)
            for iprob = 1:length(unq_subjects_probe_all)-2
                if iprob < length(unq_subjects_probe_all)-2
                    plot(diff_waits_past(isubj, [iprob iprob+1]), first_mtx(isubj, [iprob iprob+1]), '-', 'color', colors(iprob,:), 'linewidth', 0.5)
                end
                plot(diff_waits_past(isubj, iprob), first_mtx(isubj, iprob), 'o', 'color', colors(iprob,:))
            end
        end
        for iprob = 1:length(unq_subjects_probe_all)-2
            if iprob < length(unq_subjects_probe_all)-2
                plot(nanmean(diff_waits_past(:, [iprob iprob+1])), nanmean(first_mtx(:, [iprob iprob+1])), '-', 'color', colors(iprob,:), 'linewidth', 6)
            end
            plot(nanmean(diff_waits_past(:, iprob)), nanmean(first_mtx(:, iprob)), '.', 'color', colors(iprob,:), 'markersize', 60)
            fit_line(diff_waits_past(:, iprob), first_mtx(:, iprob));
        end
        set(gca,'TickLength',[0, 0]); box off;
        title(['early discrim (proportion)'])
        xlabel(['prior probe discrim'])
        ylabel(['early discrim (proportion)'])
        hold on; plot(xlim, [0 0], 'k--')
        hold on; plot([0 0], ylim, 'k--')

                % mixed model
                tbl = table(probe_hold(:), fitline_hold, subj_num_mtx(:), stage_num_mtx(:), 'VariableNames',{'dv','ProbeBeh','subject','problem'});
                tbl.subject = categorical(tbl.subject);
                tbl.problem = categorical(tbl.problem);
                lme_dissagreement_rich = fitlme(tbl, model_str)
        
        
        
        
%}
% late discrim (proportion)
%{
        % mean centered correlation dot plot
        figure; hold on; 
        var1 = diff_waits_past;
        var2 = last_mtx;
        for iprob = 1:length(unq_subjects_probe_all)-2
            % mean center
            %
            var1(:, iprob) = var1(:, iprob) - nanmean(var1(:, iprob));
            var2(:, iprob) = var2(:, iprob) - nanmean(var2(:, iprob));
            %

            %[p,r] = corr(var1(:, iprob), var2(:, iprob))

            % old
            plot(var1(:, iprob), var2(:, iprob), 'o', 'color', colors(iprob,:))
            %plot(nanmean(var1(:, iprob)), nanmean(var2(:, iprob)), '.', 'color', colors(iprob,:), 'markersize', 50)
        end
        fitline_hold_mc = var1(:);
        probe_hold_mc = var2(:,1:end-1);
        fitline_hold_mc = fitline_hold_mc(1:length(probe_hold_mc(:)));
        fitline_hold = last_mtx(:);
        probe_hold = diff_waits_past(:,1:end-1);
        fitline_hold = fitline_hold(1:length(probe_hold(:)));
        [r, p] = fit_line(fitline_hold_mc, probe_hold_mc(:));
        set(gca,'TickLength',[0, 0]); box off;
        title(['late discrim (proportion); r=' num2str(r) ', p=' num2str(p)])
        xlabel(['prior probe discrim'])
        ylabel(['late discrim (proportion)'])
        hold on; plot(xlim, [0 0], 'k--')
        hold on; plot([0 0], ylim, 'k--')
        
        
        % change over problem set dot plot
        figure; hold on;
        for isubj = 1:size(diff_waits_past,1)
            for iprob = 1:length(unq_subjects_probe_all)-2
                if iprob < length(unq_subjects_probe_all)-2
                    plot(diff_waits_past(isubj, [iprob iprob+1]), last_mtx(isubj, [iprob iprob+1]), '-', 'color', colors(iprob,:), 'linewidth', 0.5)
                end
                plot(diff_waits_past(isubj, iprob), last_mtx(isubj, iprob), 'o', 'color', colors(iprob,:))
            end
        end
        for iprob = 1:length(unq_subjects_probe_all)-2
            if iprob < length(unq_subjects_probe_all)-2
                plot(nanmean(diff_waits_past(:, [iprob iprob+1])), nanmean(last_mtx(:, [iprob iprob+1])), '-', 'color', colors(iprob,:), 'linewidth', 6)
            end
            plot(nanmean(diff_waits_past(:, iprob)), nanmean(last_mtx(:, iprob)), '.', 'color', colors(iprob,:), 'markersize', 60)
            fit_line(diff_waits_past(:, iprob), last_mtx(:, iprob));
        end
        set(gca,'TickLength',[0, 0]); box off;
        title(['late discrim (proportion)'])
        xlabel(['prior probe discrim'])
        ylabel(['late discrim (proportion)'])
        hold on; plot(xlim, [0 0], 'k--')
        hold on; plot([0 0], ylim, 'k--')

                % mixed model
                tbl = table(probe_hold(:), fitline_hold, subj_num_mtx(:), stage_num_mtx(:), 'VariableNames',{'dv','ProbeBeh','subject','problem'});
                tbl.subject = categorical(tbl.subject);
                tbl.problem = categorical(tbl.problem);
                lme_dissagreement_rich = fitlme(tbl, model_str)
                
                
%}            
% early discrim (first / last)
%

        % mean centered correlation dot plot
        figure; hold on; 
        var1 = diff_waits_past;
        var2 = first_day_d;
        for iprob = 1:length(unq_subjects_probe_all)-2
            % mean center
            %
            var1(:, iprob) = var1(:, iprob) - nanmean(var1(:, iprob));
            var2(:, iprob) = var2(:, iprob) - nanmean(var2(:, iprob));
            %}

            %[p,r] = corr(var1(:, iprob), var2(:, iprob))

            % old
            plot(var1(:, iprob), var2(:, iprob), 'o', 'color', colors(iprob,:))
            %plot(nanmean(var1(:, iprob)), nanmean(var2(:, iprob)), '.', 'color', colors(iprob,:), 'markersize', 50)
        end
        fitline_hold_mc = var1(:);
        probe_hold_mc = var2(:,1:end-1);
        fitline_hold_mc = fitline_hold_mc(1:length(probe_hold_mc(:)));
        fitline_hold = first_day_d(:);
        probe_hold = diff_waits_past(:,1:end-1);
        fitline_hold = fitline_hold(1:length(probe_hold(:)));
        [r, p] = fit_line(fitline_hold_mc, probe_hold_mc(:));
        set(gca,'TickLength',[0, 0]); box off;
        title(['early discrim (first / last); r=' num2str(r) ', p=' num2str(p)])
        xlabel(['prior probe discrim'])
        ylabel(['early discrim (first / last)'])
        hold on; plot(xlim, [0 0], 'k--')
        hold on; plot([0 0], ylim, 'k--')
        
        
        % change over problem set dot plot
        figure; hold on;
        for isubj = 1:size(diff_waits_past,1)
            for iprob = 1:length(unq_subjects_probe_all)-2
                if iprob < length(unq_subjects_probe_all)-2
                    plot(diff_waits_past(isubj, [iprob iprob+1]), first_day_d(isubj, [iprob iprob+1]), '-', 'color', colors(iprob,:), 'linewidth', 0.5)
                end
                plot(diff_waits_past(isubj, iprob), first_day_d(isubj, iprob), 'o', 'color', colors(iprob,:))
            end
        end
        
        r_and_p = nan(length(unq_subjects_probe_all)-2, 2);%all r vals and p vals
        for iprob = 1:length(unq_subjects_probe_all)-2
            if iprob < length(unq_subjects_probe_all)-2
                plot(nanmean(diff_waits_past(:, [iprob iprob+1])), nanmean(first_day_d(:, [iprob iprob+1])), '-', 'color', colors(iprob,:), 'linewidth', 6)
            end
            plot(nanmean(diff_waits_past(:, iprob)), nanmean(first_day_d(:, iprob)), '.', 'color', colors(iprob,:), 'markersize', 60)
            [r_and_p(iprob,1), r_and_p(iprob,2)] = fit_line(diff_waits_past(:, iprob), first_day_d(:, iprob));
        end
        set(gca,'TickLength',[0, 0]); box off;
        title(['early discrim (first / last)'])
        xlabel(['prior probe discrim'])
        ylabel(['early discrim (first / last)'])
        hold on; plot(xlim, [0 0], 'k--')
        hold on; plot([0 0], ylim, 'k--')

                % mixed model
                tbl = table(probe_hold(:), fitline_hold, subj_num_mtx(:), stage_num_mtx(:), 'VariableNames',{'dv','ProbeBeh','subject','problem'});
                tbl.subject = categorical(tbl.subject);
                tbl.problem = categorical(tbl.problem);
                lme_early_discrim = fitlme(tbl, model_str)
                
         % r val bar     
         figure; hold on
         bar(r_and_p(:,1)); 
         sig_asterisks(r_and_p(:,2), 1:length(r_and_p(:,2)), r_and_p(:,1).*1.2)
         ylim([-1.09 1.09])
         title(['early discrim (first / last)'])
        
        
        
%
% late discrim (first / last)
%
        % mean centered correlation dot plot
        figure; hold on; 
        var1 = diff_waits_past;
        var2 = last_day_d;
        for iprob = 1:length(unq_subjects_probe_all)-2
            % mean center
            %
            var1(:, iprob) = var1(:, iprob) - nanmean(var1(:, iprob));
            var2(:, iprob) = var2(:, iprob) - nanmean(var2(:, iprob));
            %}

            %[p,r] = corr(var1(:, iprob), var2(:, iprob))

            % old
            plot(var1(:, iprob), var2(:, iprob), 'o', 'color', colors(iprob,:))
            %plot(nanmean(var1(:, iprob)), nanmean(var2(:, iprob)), '.', 'color', colors(iprob,:), 'markersize', 50)
        end
        fitline_hold_mc = var1(:);
        probe_hold_mc = var2(:,1:end-1);
        fitline_hold_mc = fitline_hold_mc(1:length(probe_hold_mc(:)));
        fitline_hold = last_day_d(:);
        probe_hold = diff_waits_past(:,1:end-1);
        fitline_hold = fitline_hold(1:length(probe_hold(:)));
        [r, p] = fit_line(fitline_hold_mc, probe_hold_mc(:));
        set(gca,'TickLength',[0, 0]); box off;
        title(['late discrim (first / last); r=' num2str(r) ', p=' num2str(p)])
        xlabel(['prior probe discrim'])
        ylabel(['late discrim (first / last)'])
        hold on; plot(xlim, [0 0], 'k--')
        hold on; plot([0 0], ylim, 'k--')
        
        
        % change over problem set dot plot
        figure; hold on;
        for isubj = 1:size(diff_waits_past,1)
            for iprob = 1:length(unq_subjects_probe_all)-2
                if iprob < length(unq_subjects_probe_all)-2
                    plot(diff_waits_past(isubj, [iprob iprob+1]), last_day_d(isubj, [iprob iprob+1]), '-', 'color', colors(iprob,:), 'linewidth', 0.5)
                end
                plot(diff_waits_past(isubj, iprob), last_day_d(isubj, iprob), 'o', 'color', colors(iprob,:))
            end
        end
        r_and_p = nan(length(unq_subjects_probe_all)-2, 2);%all r vals and p vals
        for iprob = 1:length(unq_subjects_probe_all)-2
            if iprob < length(unq_subjects_probe_all)-2
                plot(nanmean(diff_waits_past(:, [iprob iprob+1])), nanmean(last_day_d(:, [iprob iprob+1])), '-', 'color', colors(iprob,:), 'linewidth', 6)
            end
            plot(nanmean(diff_waits_past(:, iprob)), nanmean(last_day_d(:, iprob)), '.', 'color', colors(iprob,:), 'markersize', 60)
            [r_and_p(iprob,1), r_and_p(iprob,2)] = fit_line(diff_waits_past(:, iprob), last_day_d(:, iprob));
        end
        set(gca,'TickLength',[0, 0]); box off;
        title(['late discrim (first / last)'])
        xlabel(['prior probe discrim'])
        ylabel(['late discrim (first / last)'])
        hold on; plot(xlim, [0 0], 'k--')
        hold on; plot([0 0], ylim, 'k--')

                % mixed model
                tbl = table(probe_hold(:), fitline_hold, subj_num_mtx(:), stage_num_mtx(:), 'VariableNames',{'dv','ProbeBeh','subject','problem'});
                tbl.subject = categorical(tbl.subject);
                tbl.problem = categorical(tbl.problem);
                lme_late_discrim = fitlme(tbl, model_str)
                
                
         % r val bar     
         figure; hold on
         bar(r_and_p(:,1)); 
         sig_asterisks(r_and_p(:,2), 1:length(r_and_p(:,2)), r_and_p(:,1).*1.2)
         ylim([-1.09 1.09])
         title(['late discrim (first / last)'])

%{
% learning rate (late discrim - early discrim) / number of trials
figure; hold on; 
for iprob = 1:6
    plot(diff_waits_past(:, iprob), learn_rate_mtx(:, iprob), 'o', 'color', colors(iprob,:))
    plot(nanmean(diff_waits_past(:, iprob)), nanmean(learn_rate_mtx(:, iprob)), '.', 'color', colors(iprob,:), 'markersize', 50)
end
[r, p] = fit_line(diff_waits_past(:), learn_rate_mtx(:));
set(gca,'TickLength',[0, 0]); box off;
title(['learning rate; r=' num2str(r) ', p=' num2str(p)])
xlabel(['prior probe discrim'])
ylabel(['learning rate'])
hold on; plot(xlim, [0 0], 'k--')
hold on; plot([0 0], ylim, 'k--')

        % mixed model
        tbl = table(learn_rate_mtx(:), diff_waits_past(:), subj_num_mtx(:), stage_num_mtx(:),'VariableNames',{'dv','ProbeBeh','subject','problem'});
        tbl.subject = categorical(tbl.subject);
        tbl.problem = categorical(tbl.problem);
        lme_LearnRate_past = fitlme(tbl,model_str)
%}

%{
% amount of learning (late discrim - early discrim)
figure; hold on; 
for iprob = 1:6
    plot(diff_waits_past(:, iprob), delta_mtx(:, iprob), 'o', 'color', colors(iprob,:))
    plot(nanmean(diff_waits_past(:, iprob)), nanmean(delta_mtx(:, iprob)), '.', 'color', colors(iprob,:), 'markersize', 50)
end
[r, p] = fit_line(diff_waits_past(:), delta_mtx(:));
set(gca,'TickLength',[0, 0]); box off;
title(['amount learned; r=' num2str(r) ', p=' num2str(p)])
xlabel(['prior probe discrim'])
ylabel(['amount learned'])
hold on; plot(xlim, [0 0], 'k--')
hold on; plot([0 0], ylim, 'k--')

        % mixed model
        tbl = table(delta_mtx(:), diff_waits_past(:), subj_num_mtx(:), stage_num_mtx(:),'VariableNames',{'dv','ProbeBeh','subject','problem'});
        tbl.subject = categorical(tbl.subject);
        tbl.problem = categorical(tbl.problem);
        lme_AmtLearn_past = fitlme(tbl,model_str)
%}

%{
% number of trials to some minimum discrimination
figure; hold on; 
for iprob = 1:6
    plot(diff_waits_past(:, iprob), min_discrim_mtx(:, iprob), 'o', 'color', colors(iprob,:))
    plot(nanmean(diff_waits_past(:, iprob)), nanmean(min_discrim_mtx(:, iprob)), '.', 'color', colors(iprob,:), 'markersize', 50)
end
[r, p] = fit_line(diff_waits_past(:), min_discrim_mtx(:));
set(gca,'TickLength',[0, 0]); box off;
title(['trials to min discrim; r=' num2str(r) ', p=' num2str(p)])
xlabel(['prior probe discrim'])
ylabel(['trials to min discrim(0.3)'])
hold on; plot(xlim, [0 0], 'k--')
hold on; plot([0 0], ylim, 'k--')

        % mixed model
        tbl = table(min_discrim_mtx(:), diff_waits_past(:), subj_num_mtx(:), stage_num_mtx(:),'VariableNames',{'dv','ProbeBeh','subject','problem'});
        tbl.subject = categorical(tbl.subject);
        tbl.problem = categorical(tbl.problem);
        lme_MinDiscr_past = fitlme(tbl,model_str)
%}

%% FUTURE probe plots NOTHING SIGNIFICANT
%{
figure; hold on; 
for iprob = 1:6
    plot(diff_waits_future(:, iprob), all_sesh_ct(:, iprob), 'o', 'color', colors(iprob,:))
    plot(nanmean(diff_waits_future(:, iprob)), nanmean(all_sesh_ct(:, iprob)), '.', 'color', colors(iprob,:), 'markersize', 50)
end
[r, p] = fit_line(diff_waits_future(:), all_sesh_ct(:));
set(gca,'TickLength',[0, 0]); box off;
title(['num trials to crit; r=' num2str(r) ', p=' num2str(p)])
xlabel(['next probe discrim'])
ylabel(['num trials to crit'])
hold on; plot(xlim, [0 0], 'k--')
hold on; plot([0 0], ylim, 'k--')

        % mixed model
        tbl = table(all_sesh_ct(:), diff_waits_future(:), subj_num_mtx(:), stage_num_mtx(:),'VariableNames',{'dv','ProbeBeh','subject','problem'});
        tbl.subject = categorical(tbl.subject);
        tbl.problem = categorical(tbl.problem);
        lme_NumTrials_future = fitlme(tbl,model_str)


figure; hold on; 
for iprob = 1:6
    plot(diff_waits_future(:, iprob), first_mtx(:, iprob), 'o', 'color', colors(iprob,:))
    plot(nanmean(diff_waits_future(:, iprob)), nanmean(first_mtx(:, iprob)), '.', 'color', colors(iprob,:), 'markersize', 50)
end
[r, p] = fit_line(diff_waits_future(:), first_mtx(:));
set(gca,'TickLength',[0, 0]); box off;
title(['early discrim; r=' num2str(r) ', p=' num2str(p)])
xlabel(['next probe discrim'])
ylabel(['early discrim'])
hold on; plot(xlim, [0 0], 'k--')
hold on; plot([0 0], ylim, 'k--')

        % mixed model
        tbl = table(first_mtx(:), diff_waits_future(:), subj_num_mtx(:), stage_num_mtx(:),'VariableNames',{'dv','ProbeBeh','subject','problem'});
        tbl.subject = categorical(tbl.subject);
        tbl.problem = categorical(tbl.problem);
        lme_EarlyDiscr_future = fitlme(tbl,model_str)


figure; hold on; 
for iprob = 1:6
    plot(diff_waits_future(:, iprob), last_mtx(:, iprob), 'o', 'color', colors(iprob,:))
    plot(nanmean(diff_waits_future(:, iprob)), nanmean(last_mtx(:, iprob)), '.', 'color', colors(iprob,:), 'markersize', 50)
end
[r, p] = fit_line(diff_waits_future(:), last_mtx(:));
set(gca,'TickLength',[0, 0]); box off;
title(['late discrim; r=' num2str(r) ', p=' num2str(p)])
xlabel(['next probe discrim'])
ylabel(['late discrim'])
hold on; plot(xlim, [0 0], 'k--')
hold on; plot([0 0], ylim, 'k--')

        % mixed model
        tbl = table(last_mtx(:), diff_waits_future(:), subj_num_mtx(:), stage_num_mtx(:),'VariableNames',{'dv','ProbeBeh','subject','problem'});
        tbl.subject = categorical(tbl.subject);
        tbl.problem = categorical(tbl.problem);
        lme_LateDiscr_future = fitlme(tbl,model_str)


figure; hold on; 
for iprob = 1:6
    plot(diff_waits_future(:, iprob), learn_rate_mtx(:, iprob), 'o', 'color', colors(iprob,:))
    plot(nanmean(diff_waits_future(:, iprob)), nanmean(learn_rate_mtx(:, iprob)), '.', 'color', colors(iprob,:), 'markersize', 50)
end
[r, p] = fit_line(diff_waits_future(:), learn_rate_mtx(:));
set(gca,'TickLength',[0, 0]); box off;
title(['learning rate; r=' num2str(r) ', p=' num2str(p)])
xlabel(['next probe discrim'])
ylabel(['learning rate'])
hold on; plot(xlim, [0 0], 'k--')
hold on; plot([0 0], ylim, 'k--')

        % mixed model
        tbl = table(learn_rate_mtx(:), diff_waits_future(:), subj_num_mtx(:), stage_num_mtx(:),'VariableNames',{'dv','ProbeBeh','subject','problem'});
        tbl.subject = categorical(tbl.subject);
        tbl.problem = categorical(tbl.problem);
        lme_LearnRate_future = fitlme(tbl,model_str)
        

figure; hold on; 
for iprob = 1:6
    plot(diff_waits_future(:, iprob), delta_mtx(:, iprob), 'o', 'color', colors(iprob,:))
    plot(nanmean(diff_waits_future(:, iprob)), nanmean(delta_mtx(:, iprob)), '.', 'color', colors(iprob,:), 'markersize', 50)
end
[r, p] = fit_line(diff_waits_future(:), delta_mtx(:));
set(gca,'TickLength',[0, 0]); box off;
title(['amount learned; r=' num2str(r) ', p=' num2str(p)])
xlabel(['next probe discrim'])
ylabel(['amount learned'])
hold on; plot(xlim, [0 0], 'k--')
hold on; plot([0 0], ylim, 'k--')

        % mixed model
        tbl = table(delta_mtx(:), diff_waits_future(:), subj_num_mtx(:), stage_num_mtx(:),'VariableNames',{'dv','ProbeBeh','subject','problem'});
        tbl.subject = categorical(tbl.subject);
        tbl.problem = categorical(tbl.problem);
        lme_AmtLearn_future = fitlme(tbl,model_str)


figure; hold on; 
for iprob = 1:6
    plot(diff_waits_future(:, iprob), min_discrim_mtx(:, iprob), 'o', 'color', colors(iprob,:))
    plot(nanmean(diff_waits_future(:, iprob)), nanmean(min_discrim_mtx(:, iprob)), '.', 'color', colors(iprob,:), 'markersize', 50)
end
[r, p] = fit_line(diff_waits_future(:), min_discrim_mtx(:));
set(gca,'TickLength',[0, 0]); box off;
title(['trials to min discrim; r=' num2str(r) ', p=' num2str(p)])
xlabel(['next probe discrim'])
ylabel(['trials to min discrim (0.3)'])
hold on; plot(xlim, [0 0], 'k--')
hold on; plot([0 0], ylim, 'k--')

        % mixed model
        tbl = table(min_discrim_mtx(:), diff_waits_future(:), subj_num_mtx(:), stage_num_mtx(:),'VariableNames',{'dv','ProbeBeh','subject','problem'});
        tbl.subject = categorical(tbl.subject);
        tbl.problem = categorical(tbl.problem);
        lme_MinDiscr_future = fitlme(tbl,model_str)
%}

%% probe delta tests NOTHING SIGNIFICANT FOR SOME REASON
% THIS WAS REPLACED WITH 'beyond rich' and 'beyond poor' analyses above
%{
figure; hold on; 
for iprob = 1:6
    plot(all_sesh_ct(:, iprob), probe_delta(:, iprob), 'o', 'color', colors(iprob,:))
    plot(nanmean(all_sesh_ct(:, iprob)), nanmean(probe_delta(:, iprob)), '.', 'color', colors(iprob,:), 'markersize', 50)
end
[r, p] = fit_line(all_sesh_ct(:), probe_delta(:));
set(gca,'TickLength',[0, 0]); box off;
title(['Probe delta vs num trials; r=' num2str(r) ', p=' num2str(p)])
xlabel(['num trials to crit'])
ylabel(['Probe response delta'])
%hold on; plot(xlim, [0 0], 'k--')
%hold on; plot([0 0], ylim, 'k--')

        % mixed model
        tbl = table(probe_delta(:), all_sesh_ct(:), subj_num_mtx(:), stage_num_mtx(:),'VariableNames',{'dv','ProbeBeh','subject','problem'});
        tbl.subject = categorical(tbl.subject);
        tbl.problem = categorical(tbl.problem);
        lme_numtrials = fitlme(tbl,model_str)
        
        
figure; hold on; 
for iprob = 1:6
    plot(first_mtx(:, iprob), probe_delta(:, iprob), 'o', 'color', colors(iprob,:))
    plot(nanmean(first_mtx(:, iprob)), nanmean(probe_delta(:, iprob)), '.', 'color', colors(iprob,:), 'markersize', 50)
end
[r, p] = fit_line(first_mtx(:), probe_delta(:));
set(gca,'TickLength',[0, 0]); box off;
title(['Probe delta vs early discrim; r=' num2str(r) ', p=' num2str(p)])
xlabel(['early discrim'])
ylabel(['Probe response delta'])
%hold on; plot(xlim, [0 0], 'k--')
%hold on; plot([0 0], ylim, 'k--')

        % mixed model
        tbl = table(probe_delta(:), first_mtx(:), subj_num_mtx(:), stage_num_mtx(:),'VariableNames',{'dv','ProbeBeh','subject','problem'});
        tbl.subject = categorical(tbl.subject);
        tbl.problem = categorical(tbl.problem);
        lme_firstmtx = fitlme(tbl,model_str)
        
        
figure; hold on; 
for iprob = 1:6
    plot(last_mtx(:, iprob), probe_delta(:, iprob), 'o', 'color', colors(iprob,:))
    plot(nanmean(last_mtx(:, iprob)), nanmean(probe_delta(:, iprob)), '.', 'color', colors(iprob,:), 'markersize', 50)
end
[r, p] = fit_line(last_mtx(:), probe_delta(:));
set(gca,'TickLength',[0, 0]); box off;
title(['Probe delta vs late discrim; r=' num2str(r) ', p=' num2str(p)])
xlabel(['late discrim'])
ylabel(['Probe response delta'])
%hold on; plot(xlim, [0 0], 'k--')
%hold on; plot([0 0], ylim, 'k--')

        % mixed model
        tbl = table(probe_delta(:), last_mtx(:), subj_num_mtx(:), stage_num_mtx(:),'VariableNames',{'dv','ProbeBeh','subject','problem'});
        tbl.subject = categorical(tbl.subject);
        tbl.problem = categorical(tbl.problem);
        lme_lastmtx = fitlme(tbl,model_str)
        
        
figure; hold on; 
for iprob = 1:6
    plot(learn_rate_mtx(:, iprob), probe_delta(:, iprob), 'o', 'color', colors(iprob,:))
    plot(nanmean(learn_rate_mtx(:, iprob)), nanmean(probe_delta(:, iprob)), '.', 'color', colors(iprob,:), 'markersize', 50)
end
[r, p] = fit_line(learn_rate_mtx(:), probe_delta(:));
set(gca,'TickLength',[0, 0]); box off;
title(['Probe delta vs learn rate; r=' num2str(r) ', p=' num2str(p)])
xlabel(['learn rate'])
ylabel(['Probe response delta'])
%hold on; plot(xlim, [0 0], 'k--')
%hold on; plot([0 0], ylim, 'k--')

        % mixed model
        tbl = table(probe_delta(:), learn_rate_mtx(:), subj_num_mtx(:), stage_num_mtx(:),'VariableNames',{'dv','ProbeBeh','subject','problem'});
        tbl.subject = categorical(tbl.subject);
        tbl.problem = categorical(tbl.problem);
        lme_learnrate = fitlme(tbl,model_str)
        
        
figure; hold on; 
for iprob = 1:6
    plot(delta_mtx(:, iprob), probe_delta(:, iprob), 'o', 'color', colors(iprob,:))
    plot(nanmean(delta_mtx(:, iprob)), nanmean(probe_delta(:, iprob)), '.', 'color', colors(iprob,:), 'markersize', 50)
end
[r, p] = fit_line(delta_mtx(:), probe_delta(:));
set(gca,'TickLength',[0, 0]); box off;
title(['Probe delta vs Problem delta; r=' num2str(r) ', p=' num2str(p)])
xlabel(['Problem response delta'])
ylabel(['Probe response delta'])
%hold on; plot(xlim, [0 0], 'k--')
%hold on; plot([0 0], ylim, 'k--')

        % mixed model
        tbl = table(probe_delta(:), delta_mtx(:), subj_num_mtx(:), stage_num_mtx(:),'VariableNames',{'dv','ProbeBeh','subject','problem'});
        tbl.subject = categorical(tbl.subject);
        tbl.problem = categorical(tbl.problem);
        lme_deltas = fitlme(tbl,model_str)
        
        
figure; hold on; 
for iprob = 1:6
    plot(min_discrim_mtx(:, iprob), probe_delta(:, iprob), 'o', 'color', colors(iprob,:))
    plot(nanmean(min_discrim_mtx(:, iprob)), nanmean(probe_delta(:, iprob)), '.', 'color', colors(iprob,:), 'markersize', 50)
end
[r, p] = fit_line(min_discrim_mtx(:), probe_delta(:));
set(gca,'TickLength',[0, 0]); box off;
title(['Probe delta vs Problem min discrim; r=' num2str(r) ', p=' num2str(p)])
xlabel(['min discrim'])
ylabel(['Probe response delta'])
%hold on; plot(xlim, [0 0], 'k--')
%hold on; plot([0 0], ylim, 'k--')

        % mixed model
        tbl = table(probe_delta(:), min_discrim_mtx(:), subj_num_mtx(:), stage_num_mtx(:),'VariableNames',{'dv','ProbeBeh','subject','problem'});
        tbl.subject = categorical(tbl.subject);
        tbl.problem = categorical(tbl.problem);
        lme_minDiscrim = fitlme(tbl,model_str)

%}

%% SPECIAL TESTS
%{
% mixed model
tbl = table(all_sesh_ct(:), first_mtx(:), diff_waits_past(:), subj_num_mtx(:), stage_num_mtx(:),'VariableNames',{'numtrials', 'earlyDiscrim', 'ProbeBeh', 'subject', 'problem'});
tbl.subject = categorical(tbl.subject);
tbl.problem = categorical(tbl.problem);
lme_NumTrials_combo = fitlme(tbl,'numtrials~ProbeBeh+earlyDiscrim+(1|subject)+(1|subject)+(1|problem)')


tbl = table(last_mtx(:), first_mtx(:), diff_waits_past(:), subj_num_mtx(:), stage_num_mtx(:),'VariableNames',{'lateDiscrim', 'earlyDiscrim', 'ProbeBeh', 'subject', 'problem'});
tbl.subject = categorical(tbl.subject);
tbl.problem = categorical(tbl.problem);
lme_LateDiscr_combo = fitlme(tbl,'lateDiscrim~ProbeBeh+earlyDiscrim+(1|subject)+(1|subject)+(1|problem)')
end
%}
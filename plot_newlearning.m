function plot_newlearning(datafolder, nlearn_num, day_delay, learn_days)
% plots curve using subject means at each tone frequency
% probe_num is 1 2 and/or 3 (first, second, or third probe session)
% delay day is 1 and/or 30 (delay between last training day and first
% probe)

%figure
figure; hold on

%input checks
if size(learn_days,1)>size(learn_days,2)
    learn_days = learn_days';
end
if size(nlearn_num,1)>size(nlearn_num,2)
    nlearn_num = nlearn_num';
end

%colors
%colors = distinguishable_colors(length(learn_days));
plot_lobound = min(learn_days);
plot_hibound = min([max(learn_days) 14]);
colors = parula(plot_hibound-plot_lobound+1);

%legend trick
for ilgd = 1:size(colors,1)
    plot(6500,65,'o', 'color', colors(ilgd,:))
    legend_input{ilgd} = ['d' num2str(ilgd)];
end

% folderpath
folderpath = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\';
folderpath = [folderpath datafolder];

%get new learning files
nlearn_files = get_file_paths_targeted_II(folderpath, {'newlearn'}, {[num2str(day_delay) 'd']});
if isempty(nlearn_files)
    return
end

%nlearn_files = nlearn_files(contains(nlearn_files,'634891') | contains(nlearn_files,'634912'));
%nlearn_files = nlearn_files(contains(nlearn_files,'634890') | contains(nlearn_files,'634913'))



%organize files by learning number
learning_num_cell = cell(1,length(nlearn_num));
iln_ct = 0;
for iln = nlearn_num    
    iln_ct = iln_ct+1;
    learn_num_str = ['0' num2str(iln)];
    learning_num_cell{iln_ct} = nlearn_files(contains(nlearn_files, ['n' learn_num_str]));
end

% remove sessions with too few probe trials
min_rich_probe_trials = 2; % minimum samples

%iterate through learning numbers and sessions, plotting lines
jittered_freqs = nan(length(nlearn_num), length(learn_days), 3);
iln_ct = 0;
for iln = nlearn_num
    freqs_all = [];
    iln_ct = iln_ct+1;
   for isesh = learn_days
       sesh_num = learn_days(isesh);
       sesh_num_str = num2str(sesh_num);
       if length(sesh_num_str)<2
           sesh_num_str = ['0' sesh_num_str];
       end      
       
       plot_files = learning_num_cell{iln_ct}(contains(learning_num_cell{iln_ct}, ['-' sesh_num_str]));
       
        if isempty(plot_files)
            continue
        end
        file_waits = nan(size(plot_files,2),3);
        for ifile = 1:size(plot_files,1)
            load(plot_files{ifile}); 
            
            % wait times
            mean_waits = wait_times_prep(trl_mtx,1);
            [~, aw_freq] = wait_times_prep(trl_mtx,2);
            [prob_dist, pd_freq] = rwd_prob_by_freq(medass_cell);
            rich_tones = pd_freq(prob_dist>.5);
            freqs = unique(aw_freq);
            freqs_all = unique([freqs_all; freqs]);

            if sum(ismember(aw_freq,rich_tones))<min_rich_probe_trials
                continue
            end
            
            file_waits(ifile,:) = mean_waits;
            
        end
        
        jittered_freqs(iln_ct, isesh, :) = jitter_xpos(freqs_all, nanmean(file_waits,1), freqs_all.*0.35);
        
        jf_hold = jittered_freqs(iln_ct, isesh, :); jf_hold = jf_hold(:);
        
        plot(jf_hold, nanmean(file_waits,1), '-', 'color', 0.95.*[1 1 1])
                
   end
end

%plot prediction normal curves
%try
    plot_probes_means(datafolder, [1], day_delay, nlearn_num, 1);
%catch
%end
ylim([0 50])

%iterate through learning numbers and sessions, plotting colored errorbars
all_waits_cell = cell(length(learn_days), length(nlearn_num));
all_waits_ids = cell(length(learn_days), length(nlearn_num));
sesh_z = cell(length(learn_days), length(nlearn_num));
count_iln = 0;
for iln = nlearn_num
    count_iln = count_iln + 1;
    isesh_ct= 0;
   for isesh = learn_days
       isesh_ct = isesh_ct+1;
       sesh_num = learn_days(isesh_ct);
       sesh_num_str = num2str(sesh_num);
       if length(sesh_num_str)<2
           sesh_num_str = ['0' sesh_num_str];
       end 

        plot_files = learning_num_cell{count_iln}(contains(learning_num_cell{count_iln}, ['-' sesh_num_str]));  

        if isempty(plot_files)
            continue
        elseif isesh_ct == 1
            unq_subj = [];
            for isubj = 1:size(plot_files,1)
                last_slash = strfind(plot_files{isubj},'\');
                last_slash = last_slash(end);
                unq_subj = [unq_subj;plot_files{isubj}(last_slash-8:last_slash-1)];
            end
        end
        sesh_z{isesh_ct, count_iln} = nan(1,size(unq_subj,1));
        
        file_waits = nan(size(plot_files,2),3);
        for ifile = 1:size(plot_files,1)
            
            % identify subj of current file
            for isubj = 1:size(unq_subj,1)
                if contains(plot_files{ifile},unq_subj(isubj,:))
                    subj_id = isubj;
                    continue
                end
            end
            
            load(plot_files{ifile});
            
            % wait times
            mean_waits = wait_times_prep(trl_mtx,1);
            [all_waits, aw_freq] = wait_times_prep(trl_mtx,2);
            [prob_dist, pd_freq] = rwd_prob_by_freq(medass_cell);
            rich_tones = pd_freq(prob_dist>.5);
            

            if sum(ismember(aw_freq, rich_tones))<min_rich_probe_trials
                continue
            end

            file_waits(ifile,:) = mean_waits;
            all_waits_cell{isesh_ct, count_iln} = [all_waits_cell{isesh_ct, count_iln}; file_waits(ifile,:)];
            all_waits_ids{isesh_ct, count_iln} = [all_waits_ids{isesh_ct, count_iln}; subj_id]; % load subj of current file
            sesh_z{isesh_ct, count_iln}(subj_id) = zdiff(all_waits(ismember(aw_freq, rich_tones)), all_waits(~ismember(aw_freq, rich_tones)));

        end

        

   end
   
   % find excluded days for each subject
   unq_subj = unique(cell2mat(all_waits_ids(:,count_iln)));
   subj_idx = nan(size(all_waits_ids,1), length(unq_subj));
   for isubj = 1:length(unq_subj)
      current_subj = unq_subj(isubj);
      for isesh = 1:size(all_waits_ids,1)
          if ismember(current_subj, all_waits_ids{isesh,count_iln})
            subj_idx(isesh,isubj) = 1;
          else
            subj_idx(isesh,isubj) = 0;
          end
      end
   end
   

   % rewrite all_waits_cell
   all_waits_cell_rewrite = cell(size(all_waits_cell(:,count_iln)));
   for isubj = 1:length(unq_subj)
       current_subj = unq_subj(isubj);
       sesh_idx = find(subj_idx(:,isubj)==1);

       for iidx = 1:length(sesh_idx)
           current_sesh_idx = sesh_idx(iidx);
           all_waits_cell_rewrite{iidx} = [all_waits_cell_rewrite{iidx}; all_waits_cell{current_sesh_idx,count_iln}(all_waits_ids{current_sesh_idx,count_iln}==current_subj,:)];
       end
   end
   all_waits_cell(:,count_iln) = all_waits_cell_rewrite;
   

   
   % rewrite jittered freqs
   iln_ct = 0;
   for iln = nlearn_num
    iln_ct = iln_ct+1;
       jittered_freqs_hold = nan(size(jittered_freqs(iln_ct,:,:)));
       jittered_freqs_hold(:,1:sum(~isnan(jittered_freqs(iln_ct,:,1))),:) = jittered_freqs(iln_ct,~isnan(jittered_freqs(iln_ct,:,1)),:);
       jittered_freqs(iln_ct,:,:) = jittered_freqs_hold;
   end
  
   
   
   % plot
   count_sesh = 0;
   for isesh = learn_days
       count_sesh = count_sesh+1;
       jf_hold = jittered_freqs(count_iln,count_sesh, :); jf_hold = jf_hold(:);
       
       try
            errorbar(jf_hold, nanmean(all_waits_cell{count_sesh, count_iln},1), nanstd(all_waits_cell{count_sesh, count_iln},[],1)./sqrt(sum(~isnan(all_waits_cell{count_sesh, count_iln}),1)), '.', 'markersize', 15, 'color', colors(isesh,:), 'linewidth', 1.5)
            %plot(jf_hold, nanmean(all_waits_cell{count_sesh, count_iln},1), '-', 'color', 0.95.*[1 1 1])

       catch
       end
       
   end
   
   
end

% aesthetics
set(gca,'TickLength',[0, 0]); 
set(gca, 'XScale', 'log');
box off;
ylim_hold = ylim;
ylim([0 ylim_hold(2)])
xlim([4650 38000])
xticks([5000 8500 14000 23000 35000])
xticklabels({'5000', '8500', '14000', '23000', '35000'})
xlabel('Tone Frequency (Hz)')
ylabel('Mean wait times (s)')
legend(legend_input, 'location', 'northeastoutside')


%compute session zdiff means
learn_means = nan(length(nlearn_num),length(learn_days));
learn_ses = nan(length(nlearn_num),length(learn_days));

count_iln = 0;
for iln = nlearn_num
    count_iln = count_iln + 1;

    hold_sesh_z = cell2mat(sesh_z(:,count_iln));
    for isubj = 1:size(hold_sesh_z,2)
        hold_sesh_z_prm = hold_sesh_z(~isnan(hold_sesh_z(:,isubj)),isubj);
        hold_sesh_z(:,isubj) = nan;
        hold_sesh_z(1:length(hold_sesh_z_prm),isubj) = hold_sesh_z_prm;
    end

    means = nanmean(hold_sesh_z,2);
    ses = nanstd(hold_sesh_z,[],2)./sqrt(sum(~isnan(hold_sesh_z),2));
    
   learn_means(count_iln,1:length(means)) = means;
   learn_ses(count_iln,1:length(ses)) = ses;
       
end


%plot learning curves
colors = cool(7);
figure; hold on

count=0;
for icurves = nlearn_num
    count = count+1;
    %colapse nans (removed excluded days)
    lm_hold = learn_means(count,~isnan(learn_means(count,:)));
    learn_means(count,:) = nan;
    learn_means(count,1:length(lm_hold)) = lm_hold;
    
    lses_hold = learn_ses(count,~isnan(learn_ses(count,:)));
    learn_ses(count,:) = nan;
    learn_ses(count,1:length(lses_hold)) = lses_hold;

    %plot 
    errorbar(learn_days, learn_means(count,:), learn_ses(count,:), 'color', colors(icurves,:))
    % errorbar(1:learn_days/2, nanmean(reshape(learn_means(icurves,:), 2, size(learn_means,2)/2)),...
    % nanmean(reshape(learn_ses(icurves,:), 2, size(learn_ses,2)/2)), 'color', colors(icurves,:))
end
ylim_hold = ylim;
if ylim_hold(1)>-2
    ylim_hold(1) = -2;
end
if ylim_hold(2)<4
    ylim_hold(2)=4;
end
ylim(ylim_hold)
xlim([plot_lobound-0.5 plot_hibound+0.5])
set(gca,'TickLength',[0, 0]); box off;
hold on; plot(xlim, [1 1].*0, 'k--')
legend_chars = {'NL01', 'NL02', 'NL03', 'NL04', 'NL05', 'NL06', 'NL07'};
legend(legend_chars(nlearn_num), 'location', 'northeastoutside')
xticks(learn_days); xlabel('Day')
ylabel('Preference for rich tone (zdiff)')



%first and last day only
count=0;
figure; hold on
for icurves = nlearn_num
    count = count+1;
    %colapse nans (removed excluded days)
    lm_hold = learn_means(count,~isnan(learn_means(count,:)));
    learn_means(count,:) = nan;
    learn_means(count,1:length(lm_hold)) = lm_hold;
    
    lses_hold = learn_ses(count,~isnan(learn_ses(count,:)));
    learn_ses(count,:) = nan;
    learn_ses(count,1:length(lses_hold)) = lses_hold;

    %plot 
    
    nnan_idx = ~isnan(learn_means(count,:));
    nnan_cumsum = cumsum(nnan_idx==1);
    max_nnan = min([find(nnan_cumsum==14) find(nnan_idx==1,1,'last')]);
    
    if isnan(learn_means(count,[1 max_nnan]))
        continue
    end
    
    errorbar([1 2], learn_means(count,[1 max_nnan]), learn_ses(count,[1 max_nnan]), 'color', colors(icurves,:))

end
xlim([.75 2.25])
set(gca,'TickLength',[0, 0]); box off;
xticks(1:2)
xticklabels({'First', 'Last'})
xlabel('Training day')
hold on; plot(xlim, [1 1].*0, 'k--')


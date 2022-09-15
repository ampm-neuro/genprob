function all_waits = plot_discrim_curves_I(fpaths)
%input string matching the subject file and additional string info to
%identify the desired trials (desired files should contain string).

% colors
colors = winter(length(fpaths));

% preallocate waits
all_waits = cell(length(fpaths),3);
mean_waits = nan(length(fpaths),3);
se_waits = nan(length(fpaths),3);
zdiffs = nan(length(fpaths),1);
num_trls = nan(length(fpaths),1);

% iteratively load files
unq_frq_sesh = [];
for ifile = 1:length(fpaths)
    load(fpaths{ifile}, 'trl_mtx', 'medass_cell');
    
    % number of trials
    num_trls(ifile) = size(trl_mtx,1);
    
    %zscore waits
    %trl_mtx(trl_mtx(:,3)==0,12) = zscore_mtx(trl_mtx(trl_mtx(:,3)==0,12));
    
    % compute waits
    %[mwd, wd, ~, unq_frq] = wait_times(trl_mtx, medass_cell,0);
    [mwd, ~] = wait_times_prep(trl_mtx, 1); %1 for means, 2 for all
    [wd, wd_freqs] = wait_times_prep(trl_mtx, 2); %1 for means, 2 for all
    unq_frq = unique(wd_freqs);
    unq_frq_sesh = [unq_frq_sesh; unq_frq'];
    
    % reward probability distribution
    [prob_dist, pd_freq] = rwd_prob_by_freq(medass_cell);
    rich_tones = pd_freq(prob_dist>.5);
    
    %z differences
    
    if ~isempty(wd(ismember(wd_freqs,rich_tones))) && ~isempty(wd(~ismember(wd_freqs,rich_tones)))
        zdiffs(ifile) = zdiff(wd(ismember(wd_freqs,rich_tones)), wd(~ismember(wd_freqs,rich_tones)));
    end
    % resort wait output
    %[unq_frq,frq_sort_idx] = sort(unq_frq);
    %mwd = mwd(frq_sort_idx);
    %wd = wd(frq_sort_idx);
    
    % load
    mean_waits(ifile,1:length(mwd)) = mwd;
    for ifrq = 1:length(unq_frq)
        current_frq = unq_frq(ifrq);
        all_waits{ifile,ifrq} = wd(wd_freqs==current_frq);
        se_waits(ifile,ifrq) = std(wd(wd_freqs==current_frq))./sqrt(sum(~isnan(wd(wd_freqs==current_frq))));
    end

end

% zdiff learning curve
hold on;
min_trls = 42;
gd_trl_n = find(num_trls>=min_trls);
plot(gd_trl_n, zdiffs, '-o')

hold on;
for i = 1:length(zdiffs)
    plot(gd_trl_n(i), zdiffs(i), '.', 'markersize', 20, 'color', colors(i,:))
end


plot(find(num_trls<min_trls), zdiffs(num_trls<min_trls), 'ro')
xlim([0.5 length(zdiffs)+0.5])
set(gca,'TickLength',[0, 0]); box off;
plot(xlim, [1 1].*0, 'k--')
xticks(1:length(zdiffs)); xlabel('Day')
ylabel('Preference for rich tone (zdiff)')

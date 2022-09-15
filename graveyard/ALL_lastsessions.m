function all_waits = ALL_lastsessions

% preallocate
all_waits = cell(1,41);
session_paths_last = [];

%hard code
all_unq_frq = [5000,5249,5511,5786,6074,6377,6695,7028,7379,7747,8133,8538,8964,9411,9880,10372,10890,11432,12002,12601,13229,13888,14581,15307,16070,16872,17713,18596,19523,20496,21518,22590,23716,24899,26140,27443,28811,30247,31755,33338,35000];

% path
folderpath = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\Gen_richards';

% iterate through sessions 
file_list_subjects = dir(folderpath);
file_list_subjects(1:2) = [];
for isubject = 1:length(file_list_subjects)
    current_subj = file_list_subjects(isubject).name;
    
    % session paths
    session_paths = get_file_paths([folderpath '\' current_subj]);
        
    % only last sessions
    for isesh = 1:10
        if length(isesh)<2
            isesh_str = ['train0' num2str(isesh)];
        else
            isesh_str = ['train' num2str(isesh)];
        end
        session_paths_last = [session_paths_last; session_paths(find(contains(session_paths, isesh_str)==1,1,'last'))];
    end
    
end


% iterate through last sessions
all_frqs=[];
all_waits_col=[];
for isesh = 1:length(session_paths_last)
    
   % load session
   load(session_paths_last{isesh})
   
   %load waits
   frq = unique(unq_frq, 'stable');
   [~, wd] = wait_times(trl_mtx, medass_cell,0);
   
   
   %zscore
   %{
   zmtx = [];
   for iz = 1:length(wd)
       zmtx = [zmtx; repmat(iz,size(wd{iz})) wd{iz}]; 
   end
   zmtx(:,2) = zscore_mtx(zmtx(:,2));
   for iz = 1:length(wd)
       wd{iz} = zmtx(zmtx(:,1)==iz,2);
   end
   %}
   
   for ifrq = 1:length(frq)
       %all_waits{all_unq_frq==frq(ifrq)} = [all_waits{all_unq_frq==frq(ifrq)}; wd{ifrq}];
        all_waits{all_unq_frq==frq(ifrq)} = [all_waits{all_unq_frq==frq(ifrq)}; nanmean(wd{ifrq})];
   end
end

%plot
figure; hold on
all_means = nan(size(all_unq_frq));
all_stds = nan(size(all_unq_frq));
all_lengths = nan(size(all_unq_frq));
for ifrq = 1:length(all_unq_frq)

    if ~isempty(all_waits{ifrq})
        
        if ifrq>=11 && ifrq<=19
            plot(all_unq_frq(ifrq), all_waits{ifrq}, 'o', 'color', [159, 10, 40]./255)
        else
            plot(all_unq_frq(ifrq), all_waits{ifrq}, 'o', 'color', 0.7.*[1 1 1])
        end
        
        %load to fit normal distribution
        all_waits_col = [all_waits_col; all_waits{ifrq}];
        all_frqs = [all_frqs; repmat(all_unq_frq(ifrq), size(all_waits{ifrq}))];
        
        %load errorbar
        all_means(ifrq) = nanmean(all_waits{ifrq});
        all_stds(ifrq) = nanstd(all_waits{ifrq});
        all_lengths(ifrq) = sum(~isnan(all_waits{ifrq}));
    end
end
errorbar(all_unq_frq(~isnan(all_means)), all_means(~isnan(all_means)), all_stds(~isnan(all_means))./sqrt(all_lengths(~isnan(all_means))), 'k-')

set(gca,'TickLength',[0, 0]); box off;
set(gca, 'XScale', 'log')
xlim([4500 38000])
ylim([0 45])


%plot
%{
figure; hold on
smooth_size = 3;
plot(all_unq_frq(~isnan(all_means)), smooth(all_means(~isnan(all_means)),smooth_size), 'k-', 'linewidth', 2)
plot(all_unq_frq(~isnan(all_means)), smooth(all_means(~isnan(all_means)),smooth_size) + smooth(all_stds(~isnan(all_means)),smooth_size), 'k-', 'linewidth', 1)
plot(all_unq_frq(~isnan(all_means)), smooth(all_means(~isnan(all_means)),smooth_size) - smooth(all_stds(~isnan(all_means)),smooth_size), 'k-', 'linewidth', 1)
%plot(all_unq_frq(~isnan(all_means)), all_means(~isnan(all_means)) + all_stds(~isnan(all_means))./sqrt(all_lengths(~isnan(all_means))), 'k--', 'linewidth', 1)
%plot(all_unq_frq(~isnan(all_means)), all_means(~isnan(all_means)) - all_stds(~isnan(all_means))./sqrt(all_lengths(~isnan(all_means))), 'k--', 'linewidth', 1)
set(gca,'TickLength',[0, 0]); box off;
set(gca, 'XScale', 'log')
%}

% overlay fit normal distribution
%{
[all_frqs_unq, ~, all_frqs_num] = unique(all_frqs);
[mu, sigma, multiplier, intercept] = ampm_normal_fit(all_frqs_num, all_waits_col);
fit_yvals = (normpdf(all_frqs_num, mu, sigma).*multiplier)+intercept;
line(all_frqs_unq(all_frqs_num), fit_yvals)


fit_yvals = (normpdf(all_frqs_num, mu, sigma).*multiplier)+intercept;
line(all_frqs_unq(all_frqs_num), fit_yvals)
%}





%plot
%{
all_counts = nan(60,length(all_unq_frq));
for ifrq = 1:length(all_unq_frq)   
    if ~isempty(all_waits{ifrq})
        all_counts(:,ifrq) = histcounts(all_waits{ifrq},0:1:60);
    end
end
all_counts_hm = all_counts;
%all_counts_hm(all_counts_hm==0) = nan;
%all_counts_hm = zscore_mtx(all_counts_hm);
all_counts_hm = norm_mtx(all_counts_hm);

figure; hold on
imagesc(all_counts_hm)
axis([0.5 length(all_unq_frq)+0.5 0.5 60.5])
%errorbar(1:length(all_unq_frq), inpaint_nans(all_means), all_stds./sqrt(all_lengths), 'k-')
%}






function [h1, unq_frq, mean_wait_times, coefEsts] = wait_times_plot(trl_mtx, varargin)
%plots the mean wait times at each frequency

if nargin == 2
    plot_type = varargin{1};
    current_color = [rand(1) rand(1) rand(1)];
elseif nargin == 3
    plot_type = varargin{1};
    current_color = varargin{2};
else
    plot_type = 1;
    %current_color = [rand(1) rand(1) rand(1)].*0.9;
    current_color = [0 0 0 ];
end


%colors
rich_red = [166 28 28]./255;
poor_gray = [81 81 81]./255;

% figure
hold on;

% preallocate coefEsts
coefEsts = nan(1,4);

% all wait times and corresponding frequencies
if length(unique(floor(trl_mtx(~isnan(trl_mtx(:,2)),2))))>10
    [wait_times, frequencies] = wait_times_prep(trl_mtx, 2); %probe
    if length(unique(frequencies))>2
        try
        	[wait_times, ~, freq_numbers] = wait_times_prep(trl_mtx, 2, 0);
            [~, coefEsts, modelFun] = ampm_normal_logistic_fit_algo(freq_numbers, wait_times);
        catch
        end
    end
else
    [wait_times, frequencies] = wait_times_prep(trl_mtx, 2);
    if length(unique(frequencies))>2
        [wait_times, ~, freq_numbers] = wait_times_prep(trl_mtx, 2, 0);
        [~, coefEsts, modelFun] = ampm_normal_logistic_fit_algo(freq_numbers, wait_times);
    end
end


% ttest
if length(unique(frequencies)) == 2
    [~, pval] = ttest2(wait_times(frequencies==frequencies(1)), wait_times(frequencies==frequencies(end)));
end

% all frequencies
load('unqfrq41', 'unqfrq41') 

% means and ses
unq_frq = unique(frequencies);
mean_wait_times = nan(1,length(unqfrq41));
se_wait_times = nan(size(mean_wait_times));
wt = cell(1,length(unqfrq41));
for ifrq = 1:length(unqfrq41)
    wt{ifrq} = wait_times(frequencies==unqfrq41(ifrq));
    mean_wait_times(ifrq) = mean(wait_times(frequencies==unqfrq41(ifrq)));
    se_wait_times(ifrq) = std(wait_times(frequencies==unqfrq41(ifrq)))/sqrt(sum(frequencies==unqfrq41(ifrq)));
end



% plot individual waits
if ismember(plot_type, [1 3])
     
    for ifrq = 1:length(unq_frq)
            if length(wait_times(frequencies==unq_frq(ifrq)))>1
                xpos_jit = jitter_xpos(repmat(unq_frq(ifrq), size(wait_times(frequencies==unq_frq(ifrq)))),wait_times(frequencies==unq_frq(ifrq)), unq_frq(ifrq)/2);
            else
                xpos_jit = repmat(unq_frq(ifrq), size(wait_times(frequencies==unq_frq(ifrq))));
            end
        if ifrq==1
            h1 = plot(xpos_jit, wait_times(frequencies==unq_frq(ifrq)), 'o', 'color', current_color);
        else
            plot(xpos_jit, wait_times(frequencies==unq_frq(ifrq)), 'o', 'color', h1.Color);
        end
    end
end

% plot means and se
hold on



if ismember(plot_type, [2 3])
   errorbar(unq_frq, mean_wait_times(ismember(unqfrq41, unq_frq)), se_wait_times(ismember(unqfrq41, unq_frq)), 'linewidth', 1.5, 'color', current_color);
end



% aesthetics
set(gca,'TickLength',[0, 0]);
xlabel('Tone Frequency (Hz)')
ylabel('Wait Durations (s)')

%ylim_hold = ylim; ylim([0 ylim_hold(2)]);
set(gca, 'XScale', 'log')
if length(unique(frequencies))==2
    title(['Waits; pval = ' num2str(pval)])
else
    title('Waits')
    
    hold on
    try
        plot(unqfrq41, modelFun(coefEsts, 1:41), 'k-');
    catch
    end
    
end
xlim([4500 38500])
xticks([5000 8500 14000 23000 35000])



end

    




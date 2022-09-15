function [mus, sigmas] = plot_probes_each_subj(probe_num, day_delay, varargin)
% plots curve using subject means at each tone frequency
% probe_num is 1 2 and/or 3 (first, second, or third probe session)
% delay day is 1 and/or 30 (delay between last training day and first
% probe)

%input checks
if size(probe_num,1)>size(probe_num,2)
    probe_num = probe_num';
end
if size(day_delay,1)>size(day_delay,2)
    day_delay = day_delay';
end
if nargin ==3
    plot_on = varargin{1};
elseif nargin == 4
    
else
    plot_on = 1; %default is plot ON
end

%output
mus = [];
sigmas = [];

%where files are stored
data_files_folder = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\Gen_richards';

%get probe files
all_files = get_file_paths_all(data_files_folder)';
probe_files_all = all_files(contains(all_files, 'probe'));

%unique subjects
unq_subj = [];
for ipf = 1:size(probe_files_all,2)
    start_pt = length(data_files_folder)+2;  
    end_pt = strfind(probe_files_all{ipf}(start_pt:end), '\'); 
        end_pt = start_pt + end_pt(1) - 2;
    subj_text = probe_files_all{ipf}(start_pt:end_pt);
    unq_subj{ipf} = subj_text;
end
unq_subj = unique(unq_subj)';


%probe num index
num_idx = false(size(probe_files_all));
for ipn = probe_num
    probe_num_str = num2str(ipn);
    if length(probe_num_str)<2
        probe_num_str = ['0' probe_num_str]; 
    end
    num_idx = num_idx | contains(probe_files_all, ['-' probe_num_str]);
end

%reduce
probe_files_constrained = probe_files_all(num_idx);
if isempty(probe_files_constrained)
    display('requested probe days do not exist'); return
end


%only include files from subjects that experienced the input delay day condition(s)
delay_num_subjs = [];
for iunq = 1:length(unq_subj)
    for idd = day_delay
        
        day_delay_str = num2str(idd); 
        if length(day_delay_str)<2
            day_delay_str = ['0' day_delay_str]; 
        end

        subj_idx = contains(probe_files_all, unq_subj{iunq});
        delay_num_idx = contains(probe_files_all, [day_delay_str 'd']);

        if ~isempty(probe_files_all(subj_idx & delay_num_idx))
            delay_num_subjs = [delay_num_subjs; unq_subj(iunq)];
        end
    
    end
end

%reduce
if isempty(delay_num_subjs)
    display('requested probe days do not exist'); return
end
probe_files_constrained = probe_files_constrained(contains(probe_files_constrained, delay_num_subjs));

%legend prep
colors = distinguishable_colors(size(probe_files_constrained,2));
legend_input = cell(1,size(probe_files_constrained,2));
for isesh = 1:size(probe_files_constrained,2)
    load(probe_files_constrained{isesh});
    [slash_locs] = strfind(probe_files_constrained{isesh},'\');
    legend_input{isesh} = probe_files_constrained{isesh}(slash_locs(end-1)+1 : slash_locs(end)-1);
    plot([5001 5001],[0 0],'-', 'linewidth', realmin, 'color', colors(isesh,:))
end

%plot probe files
all_waits = nan(size(probe_files_constrained,1), 41);
all_fitvals = [];
for isesh = 1:size(probe_files_constrained,2)
    load(probe_files_constrained{isesh});
    [~,wd] = wait_times(trl_mtx, medass_cell,0);
    if plot_on ==1
        [coefEsts, fit_yvals, ~, fit_error] = plot_normal_fit_subj(trl_mtx, 3, colors(isesh,:));
        mus = [mus; coefEsts(1)];
        sigmas = [sigmas; coefEsts(2)];
        all_fitvals = [all_fitvals; fit_yvals];
    end
end

%means and ses of all probe fits
mean_fit = nanmean(all_fitvals,1);
std_fit = nanstd(all_fitvals,[],1);
se_fit = std_fit./sqrt(sum(~isnan(all_fitvals)));


if plot_on ==1
    %errorbar(unique(unq_frq, 'stable'), nanmean(all_waits,1), nanstd(all_waits,[],1)./sqrt(sum(~isnan(all_waits),1)), 'k-')
    set(gca,'TickLength',[0, 0]); 
    set(gca, 'XScale', 'log');
    box off;
    ylim([0 60]);
    xlim([4650 38000])
    xticks([5000 8500 14000 23000 35000])
    xticklabels({'5000', '8500', '14000', '23000', '35000'})
    xlabel('Tone Frequency (Hz)')
    ylabel('Mean wait times (s)')
    rich_bounds; %script plotting red lines around rich area
    legend(legend_input, 'location', 'northeastoutside')
    
    
end






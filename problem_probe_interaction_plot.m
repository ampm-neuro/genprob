function problem_probe_interaction_plot(training_group, probe_nums)
% plots three versions of each probe overlayed with problem behavior
% 1. overlayed with last day of last problem
% 2. overlayed with first day of current problem
% 3. overlayed with last day of current problem


%% all probe wait times

% smooth wait times from every probe for every subject (probe,tone,subj)
[probe_wait_times, unq_subjects_probe_all] = ALL_probe_wait_times(training_group, 2);

% number of probes
num_probes = size(probe_wait_times,1);
num_tones = size(probe_wait_times,2);
num_subjects = size(probe_wait_times,3);



%% Probe wait time plots

% iterate through probes
probe_nums = flipud(sort(probe_nums));

for iprobe = probe_nums
    titlestring = ['Probe ' num2str(iprobe)];

    
    % three iterations
    for iter = 4%:-1:1
        
        % exceptions
        if iprobe >= 7 && ismember(iter, [2 3])
            continue
        end
    
        figure; hold on
        
        % wait times on this probe from all subjects (subj, tone)
        probes_local = squeeze(probe_wait_times(iprobe, :, :))';
    
        % plot average wait times at each frequency
        for isubj = 1:num_subjects 
            plot(probes_local(isubj,:),'color', [102 178 204]./255)
        end
        errorbar(nanmean(probes_local,1), nanstd(probes_local,[],1)./sqrt(sum(~isnan(probes_local),1)), 'linewidth', 2, 'color', [0 124 204]./255);
        
        % plot problem behavior
        if iter == 1 && iprobe < 8 && iprobe >1
            plot_all_subj_problem(training_group, iprobe-1, 2, 2)
            titlestring = ['Probe ' num2str(iprobe) '; last day of problem ' num2str(iprobe-1)];
        elseif iter == 2
            plot_all_subj_problem(training_group, iprobe, 1, 2)
            titlestring = ['Probe ' num2str(iprobe) '; first day of problem ' num2str(iprobe)];
        elseif iter == 3 
            plot_all_subj_problem(training_group, iprobe, 2, 2)
            titlestring = ['Probe ' num2str(iprobe) '; last day of problem ' num2str(iprobe)];
        elseif iter == 4 && iprobe ~=7 
            %close
            %continue
        end
        
        % aesthetics
        ylim([0 40])
        ylabel('Wait time (s)')
        xlim([-1 43])
        xlabel('Tone')
        set(gca,'TickLength',[0, 0]); box off;
        title(titlestring)
        xticks([1 41])
        xticklabels([5000 35000])
        
        %}
        
        % save
        print(['C:\Users\ampm1\Documents\LabMeet\LabMeet07\probe_aggregate_plot\' titlestring], '-dpdf', '-painters', '-bestfit')
    
    end
end
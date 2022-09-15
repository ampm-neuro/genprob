function plot_all_subj_problem(training_group, problem_number, first_last, freqs_or_nums)
% plots one rich and one poor dot per subject
% first = 1, last = 2
% freqs = 1, nums = 2

% folder
folderpath = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\';
folderpath = [folderpath training_group];

% gen file number
if problem_number <4
    gfn = ['gen0' num2str(problem_number +6)];
else
    gfn = ['gen' num2str(problem_number +6)];
end

% all paths
all_paths = get_file_paths_targeted(folderpath, gfn);

% unique subjects
unq_subj = [];
for ipath = 1:size(all_paths, 1)
    unq_subj = [unq_subj; find_subj_id(all_paths{ipath})];
end
unq_subj = unique(unq_subj, 'rows');

% preallocate
wait_means = nan(size(unq_subj,1),2);
freqs = nan(size(unq_subj,1),2);
freq_nums = nan(size(unq_subj,1),2);

% iterate through subjects
for isubj = 1:size(unq_subj,1)
    
    % just subject paths
    subj_paths = all_paths(contains(all_paths, unq_subj(isubj,:)));

    % if at least two sessions
    if first_last == 2 && size(subj_paths,1)>=2
        
        % last session
        load(subj_paths{end}, 'trl_mtx')

        % rich and poor means
        [wait_means(isubj,:), freqs(isubj,:), freq_nums(isubj,:)] = rich_and_poor(trl_mtx);
        
    elseif first_last == 1
        
        % first session
        load(subj_paths{1}, 'trl_mtx')
    
        % rich and poor means
        [wait_means(isubj,:), freqs(isubj,:), freq_nums(isubj,:)] = rich_and_poor(trl_mtx);
        
    else
        continue
    
    end
    
    
   
end


% prepare errorbar cell
wait_means_cell = cell(1,size(wait_means,2));
for icol = 1:length(wait_means_cell)
    wait_means_cell{icol} = wait_means(:,icol);
end

% plot
if freqs_or_nums == 1
    
    if contains(training_group, 'hivar')
        if first_last == 1
            errorbar_plot(wait_means_cell, 1, nanmean(freqs), [.8 .8 .8], [0 0 1], '--');
        else
            errorbar_plot(wait_means_cell, 1, nanmean(freqs), [.8 .8 .8], [0 0 1], '-');
        end
    else contains(training_group, 'mevar')
        if first_last == 1
            errorbar_plot(wait_means_cell, 1, nanmean(freqs), [.8 .8 .8], [0 1 0], '--');
        else
            errorbar_plot(wait_means_cell, 1, nanmean(freqs), [.8 .8 .8], [0 1 0], '-');
        end
    end
    
    axis([4500 38500 0 30]); axis square
    
elseif freqs_or_nums == 2
    
    if contains(training_group, 'hivar')
        if first_last == 1
            errorbar_plot(wait_means_cell, 1, nanmean(freq_nums), [.8 .8 .8], [0 0 1], '--');
        else
            errorbar_plot(wait_means_cell, 1, nanmean(freq_nums), [.8 .8 .8], [0 0 1], '-');
        end
    else contains(training_group, 'mevar')
        if first_last == 1
            errorbar_plot(wait_means_cell, 1, nanmean(freq_nums), [.8 .8 .8], [0 1 0], '--');
        else
            errorbar_plot(wait_means_cell, 1, nanmean(freq_nums), [.8 .8 .8], [0 1 0], '-');
        end
    end
    
    axis([0 42 0 30]); axis square
    
end



end




% internal function
function [rich_poor_means, rich_poor_freqs, rich_poor_freq_nums] = rich_and_poor(trl_mtx)
    % all 
    [wait_durations, freqs, freq_nums] = wait_times_prep(trl_mtx, 4);
    
    % identify rich tone
    if sum(freqs == freqs(1)) > floor(length(freqs)/2)
        
        % tone freq
        rich_poor_freqs(1) = freqs(end);
        rich_poor_freqs(2) = freqs(1);
        
        % tone number
        rich_poor_freq_nums(1) = freq_nums(end);
        rich_poor_freq_nums(2) = freq_nums(1);
        
    else
        
        % tone freq
        rich_poor_freqs(1) = freqs(1);
        rich_poor_freqs(2) = freqs(end);
        
        % tone number
        rich_poor_freq_nums(1) = freq_nums(1);
        rich_poor_freq_nums(2) = freq_nums(end);
    end

    % rich mean
    rich_poor_means(1) = nanmean(wait_durations(freqs==rich_poor_freqs(1)));
    
    % poor mean
    rich_poor_means(2) = nanmean(wait_durations(freqs~=rich_poor_freqs(1)));
    

end





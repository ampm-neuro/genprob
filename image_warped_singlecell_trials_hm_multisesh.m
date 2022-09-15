function all_session_fr_vects = image_warped_singlecell_trials_hm_multisesh(neuron, cell_regist_mtx, session_cell, tses)
% plots every trial from a single neuron

%% inputs

% all sessions that feature this neuron
session_numbers = find(cell_regist_mtx(neuron, :)>0);

%% compute

% figure
figure
sgtitle(['Neuron ' num2str(neuron)])

% soft preallocate all firing rate vectors for MDS plot
all_session_fr_vects = [];
all_session_fr_vects_indices = [];

% iterate through sessions
for isesh = 1:length(session_numbers)
    
    % load session
    current_session = session_numbers(isesh);
    load(session_cell{current_session})
    
    
    % firing rates
    %
    
        % ready figure subplot
        subplot(2, length(session_numbers), isesh)
        hold on

        % neuron number in this session
        neuron_sesh = cell_regist_mtx(neuron, current_session);

        % plot
        session_fr_vects = tw_activity_plot_trial_hm_warp(trl_mtx, trl_idx, frame_times, traces, neuron_sesh, 1:size(trl_mtx,1), [5 25], tses);

        % aesthetics
        caxis([0 1])
        colorbar off
        xlim([240 1760])
        title(session_title(current_session))
        yticks auto
        yticklabels auto
        set(gca, 'YDir','reverse')
        
        
    % behavior
    %
    
        % ready figure subplot
        subplot(2, length(session_numbers), length(session_numbers)+isesh)
        hold on
        
        % plot
        plot_trials_trlmtx(trl_mtx)
        
        % aesthetics
        set(gca, 'YDir','reverse')
        
        
    % load firing rate vectors for mds plot
    %
    
        % firing rates        
        all_session_fr_vects = [all_session_fr_vects; session_fr_vects];
               
        % indices
        session_idx = repmat(current_session, size(session_fr_vects,1), 1);
        all_session_fr_vects_indices = [all_session_fr_vects_indices; session_idx];
        
end

% draw fig
set(gcf, 'Position', [145 106 1700 796])
drawnow                            
                           
figure; imagesc(all_session_fr_vects)


%% Multidimensional scaling figure
%{
% remove trials without data
all_session_fr_vects_indices(sum(~isnan(all_session_fr_vects),2)==0) = [];
all_session_fr_vects(sum(~isnan(all_session_fr_vects),2)==0, :) = [];

% remove rwd time bins
all_session_fr_vects(:, sum(isnan(all_session_fr_vects),1)>0) = [];

% remove time bins after end of random delay
eod = cumsum(tses);
eod = eod(4)/eod(end);
eod = ceil(eod*size(all_session_fr_vects,2));
all_session_fr_vects(:, eod:end) = [];

% zscore
all_session_fr_vects = all_session_fr_vects-mean(all_session_fr_vects(:));
all_session_fr_vects = all_session_fr_vects./std(all_session_fr_vects(:));

% rows for each trial
all_session_fr_vects = all_session_fr_vects';

% correlation matrix
corr_mtx = correlate_columns(all_session_fr_vects);

% plot correlation matrix

figure; ampm_pcolor(1-corr_mtx); drawnow
xticks auto
xticklabels auto
yticks auto
yticklabels auto
axis square

% distance matrix
dist_mtx = nan(size(all_session_fr_vects,2));
for itrl_1 = 1:size(all_session_fr_vects,2)
    for itrl_2 = 1:size(all_session_fr_vects,2)
        trl_1 = all_session_fr_vects(:, itrl_1)';
        trl_2 = all_session_fr_vects(:, itrl_2)';
        dist_mtx(itrl_1, itrl_2) = pdist([trl_1; trl_2])./sqrt(length(trl_1));
    end
end
%}
%{
% make dissimilarity matrix
corr_mtx = 1-corr_mtx;
corr_mtx(logical(eye(size(corr_mtx,1)))) = 0;

% plot
[mds_coords, stress] = mds_plot_2(dist_mtx, 3, all_session_fr_vects_indices);

% draw fig
drawnow
%}

end

function title_str = session_title(session_number)

    % session ids
    probe_sessions = [1:3:19 20];
    first_problems = nan(1,6);
    last_problems = nan(1,6);
    for iprobe = 1:6
        first_problems(iprobe) = probe_sessions(iprobe)+1;
        last_problems(iprobe) = probe_sessions(iprobe)+2;
    end
    
    if ismember(session_number, probe_sessions)
        title_str = [num2str(session_number) ' (probe '  num2str(find(ismember(probe_sessions, session_number))) ')']; 
        
    elseif ismember(session_number, first_problems)
        title_str = [num2str(session_number) ' (problem ' num2str(find(ismember(first_problems, session_number))) ', first)']; 
            
    elseif ismember(session_number, last_problems)
        title_str = [num2str(session_number) ' (problem ' num2str(find(ismember(last_problems, session_number))) ', last)']; 

    else
        title_str = num2str(session_number);
    end

end



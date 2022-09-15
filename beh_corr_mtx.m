function [bcm, beh_corr_cells] = beh_corr_mtx(datafolder)
% creates corr mtx heatmap of two types of behavior variable determined by
% input.
%
% for each of 6 learning stages: first z, days to crit, last z (18 x 18)


%% gather session files

%could be input
folderpath = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\';
folderpath = [folderpath datafolder '\'];

% file list
session_files1 = get_file_paths_targeted(folderpath, {'mevar0'});
session_files2 = get_file_paths_targeted(folderpath, {'hivar0'});
session_files3 = get_file_paths_targeted(folderpath, {'ctl0'});
session_files = [session_files1; session_files2; session_files3];

% find unique subjects
subjs = [];
for isf = 1:size(session_files,1)
    subjs = [subjs; ...
        {session_files{isf}(strfind(session_files{isf}, '\\')+2 : strfind(session_files{isf}, '\\')+9)}];
end
subjs = unique(subjs);

%% compute session file stats

% preallocate
beh_corr_cells = cell(1,18); % first crit last repeated*6
for ic = 1:length(beh_corr_cells)
	beh_corr_cells{ic} = nan(size(subjs,1),1);
end

% indices
first_z_idx = 1:3:16;
crit_day_idx = 2:3:17;
last_z_idx = 3:3:18;

% iterate through each subj
for isubj = 1:size(subjs,1)
    
    % all subject session files
    sf_subj = session_files(contains(session_files, subjs{isubj}));
    
    % iterate through discirmination problems 1:6
    for iprob = 1:6
        
        % all subject session problem files
        sf_subj_prob = sf_subj(contains(sf_subj, ['r0' num2str(iprob)]) | contains(sf_subj, ['ctl0' num2str(iprob)]));
        
        if isempty(sf_subj_prob)
            continue
        end
        
        % compute values
        %
        
        % load frequencies
        load('unqfrq41', 'unqfrq41')
        
        
            % first z
            %
                % load session file
                load(sf_subj_prob{1}, 'trl_mtx')
                
                % wait times
                [wait_times, frequencies] = wait_times_prep(trl_mtx, 2);
                if contains(sf_subj_prob{1}, 'mevar')
                    rich_tones = unqfrq41(13:29);
                elseif contains(sf_subj_prob{1}, 'hivar')
                    rich_tones = unqfrq41([5 13 21 29 37]);
                end
                rich_waits = wait_times(ismember(frequencies,rich_tones));
                poor_waits = wait_times(~ismember(frequencies,rich_tones));
                
                % difference between means in terms of std
                zdiff_first = computeCohen_d(rich_waits, poor_waits);
                
                % load
                celln = first_z_idx(iprob);
                beh_corr_cells{celln}(isubj) = zdiff_first;
                    
            % days to crit
            %
                % iterate through sessions counting until first success
                sesh_count = 0;
                for isf_subj_prob = 1:size(sf_subj_prob,1)
                    sesh_count = sesh_count + 1;
                    
                    % load session file
                    load(sf_subj_prob{isf_subj_prob}, 'trl_mtx')

                    % wait times
                    [wait_times, frequencies] = wait_times_prep(trl_mtx, 2);
                    if contains(sf_subj_prob{1}, 'mevar')
                        rich_tones = unqfrq41(13:29);
                    elseif contains(sf_subj_prob{1}, 'hivar')
                        rich_tones = unqfrq41([5 13 21 29 37]);
                    end
                    rich_waits = wait_times(ismember(frequencies,rich_tones));
                    poor_waits = wait_times(~ismember(frequencies,rich_tones));
                    
                    % compute for crit requirements
                    [~,pval,~,stats] = ttest2(rich_waits, poor_waits);
               
                    % check for pass
                    if pval<0.01 && stats.tstat>0
                        
                        % load NEGATIVE days to crit
                        celln = crit_day_idx(iprob);
                        beh_corr_cells{celln}(isubj) = -sesh_count;
                        
                        % stop loading sessions
                        break
                    end
                end

            % last z
            %
                % load session file
                load(sf_subj_prob{end}, 'trl_mtx')
                
                % wait times
                [wait_times, frequencies] = wait_times_prep(trl_mtx, 2);
                if contains(sf_subj_prob{1}, 'mevar')
                    rich_tones = unqfrq41(13:29);
                elseif contains(sf_subj_prob{1}, 'hivar')
                    rich_tones = unqfrq41([5 13 21 29 37]);
                end
                rich_waits = wait_times(ismember(frequencies,rich_tones));
                poor_waits = wait_times(~ismember(frequencies,rich_tones));
                
                % difference between means in terms of std
                zdiff_last = computeCohen_d(rich_waits, poor_waits);
                
                % load
                celln = last_z_idx(iprob);
                beh_corr_cells{celln}(isubj) = zdiff_last;
    
    end
end


%% remove subject mean from each dv
bcc = cell2mat(beh_corr_cells);
for dv = 1:3
    bcc(:, dv:3:end) = bcc(:, dv:3:end) - nanmean(bcc(:, dv:3:end),2);
end
for i = 1:length(beh_corr_cells)
    beh_corr_cells{i} = bcc(:,i);
end



%% compute correlation matrix

%beh_corr_cells = cell(18,18); % first crit last repeated*6
bcm = nan(length(beh_corr_cells));
for ic = 1:length(beh_corr_cells)
    for ir = 1:length(beh_corr_cells)
        nnan_idx = ~isnan(beh_corr_cells{ic}) & ~isnan(beh_corr_cells{ir});
        c1 = beh_corr_cells{ic}(nnan_idx);
        c2 = beh_corr_cells{ir}(nnan_idx);

        %{
        if ic==17 & ir==16
            figure; fit_line(c1, c2)
        end
        %}
        
        if ~isempty(c1) && ~isempty(c2)
        	bcm(ic,ir) = corr(c1, c2);
        end
    end
end

% plot
figure; hold on;
imagesc((bcm))
%imagesc((bcm(1:3:end, 1:3:end)))
%imagesc((bcm(1:3:end, 1:3:end)))
%imagesc((bcm(sort([1:3:18 2:3:18]) , sort([1:3:18 2:3:18]))))
plot([0.5 18.5], [(0.5 : 3 : 18.5)' (0.5 : 3 : 18.5)'], 'k-', 'linewidth', 2)
plot([(0.5 : 3 : 18.5)' (0.5 : 3 : 18.5)'], [0.5 18.5], 'k-', 'linewidth', 2)
axis([0.5 18.5 0.5 18.5])
axis square; set(gca, 'YDir','reverse')
colorbar
caxis([-1 1])
xticks(1:18); yticks(1:18)
ticklabels = {'1-first', '1-time', '1-last',...
    '2-first', '2-time', '2-last', '3-first', '3-time', '3-last',...
    '4-first', '4-time', '4-last', '5-first', '5-time', '5-last',...
    '6-first', '6-time', '6-last'};
xticklabels(ticklabels); xtickangle(45)
yticklabels(ticklabels)
end




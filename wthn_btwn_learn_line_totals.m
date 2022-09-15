function wthn_btwn_learn_line_totals(datafolder, delay_day)
%plot within and between sesh learning cumulative over every learning stage
%plot sum total learning for each new learning stage

% plot colors
colors = cool(7);

% folderpath
folderpath = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\';
folderpath = [folderpath datafolder];



% for each new learning condition
num_training_stages = 10;
total_training = nan(num_training_stages, 2, 7); % training stage, mean/se, NL condition
total_newlearn = nan(7,2); % NL condition, mean/se
for inl = 1:7

    % load all session files
    %session_files_all = get_file_paths_targeted_II(folderpath, {'train'}, {[num2str(delay_day) 'd'], ['newlearn0' num2str(inl)]});
    session_files_all = get_file_paths_targeted(folderpath, {'train'}, {['newlearn0' num2str(inl)]});
    if isempty(session_files_all); continue; end
    
    % training stages
    %
    for istage = 1:num_training_stages

        % training stage string
        istage_string = num2str(istage);
        if length(istage_string)<2
            istage_string = ['0' istage_string];
        end

        % load all session files
        [session_files_stage, included_subjects] = get_file_paths_targeted_II(folderpath, {['train' istage_string]}, {[num2str(delay_day) 'd'], ['newlearn0' num2str(inl)]});
        
        if isempty(session_files_stage); continue; end

        % iterate through subjects
        total_learn_train = nan(size(included_subjects,1),2);
        for isubj = 1:size(included_subjects,1)

            % all subject files
            [session_files_subj] = get_file_paths_targeted_II(folderpath, {['train' istage_string], included_subjects(isubj,:)}, {[num2str(delay_day) 'd']});
            if isempty(session_files_subj)
                continue
            end

            % first and last session
            sesh_ct = 0;
            for isesh = [1 length(session_files_subj)]
                sesh_ct = sesh_ct+1;

                % session wait times
                [rich_waits, poor_waits] = rich_poor_waits(session_files_subj{isesh});
                if any(isnan([rich_waits; poor_waits])) && isesh==1 
                    try
                        [rich_waits, poor_waits] = rich_poor_waits(session_files_subj{isesh+1});
                    catch
                    end
                end

                %load session zdiff
                total_learn_train(isubj, sesh_ct) = zdiff(rich_waits, poor_waits);

            end
        end

        % delta zdiff
        total_learn_stage = total_learn_train(:,2)-total_learn_train(:,1);
        total_training(istage,:,inl) = [nanmean(total_learn_stage) nanstd(total_learn_stage)./sqrt(sum(~isnan(total_learn_stage)))];
    end
    
       
    
    % new learning condition
    %
    [session_files_NL, included_subjects_NL] = get_file_paths_targeted_II(folderpath, {'new'}, {[num2str(delay_day) 'd'], ['newlearn0' num2str(inl)]});
    if isempty(session_files_NL); continue; end

    % iterate through subjects
    total_learn_NL = nan(size(included_subjects_NL,1),2);
    for isubj_NL = 1:size(included_subjects_NL,1)

        % all subject files
        session_files_subj = session_files_NL(contains(session_files_NL,included_subjects_NL(isubj_NL,:)));
        
        % screen subject files for too few rich probes
        min_rich_probes = 2;
        sesh_files_idx = true(size(session_files_subj,1),1);
        for isesh = 1:size(session_files_subj,1)
            load(session_files_subj{isesh}, 'trl_mtx')
            rich_tone = unique(trl_mtx(:,2)); rich_tone = rich_tone(2);
            if sum(trl_mtx(:,2)==rich_tone & trl_mtx(:,3)==0) < min_rich_probes
                sesh_files_idx(isesh) = false;
            end
        end
        session_files_subj = session_files_subj(sesh_files_idx);

        % first and last session
        sesh_ct = 0;
        
        force_last_sesh = 14;
        last_sesh = size(session_files_subj,1);
        if last_sesh > force_last_sesh
            last_sesh = force_last_sesh;
        end
        
        for isesh = [1 last_sesh]
            sesh_ct = sesh_ct+1;

            % session wait times
            [rich_waits, poor_waits] = rich_poor_waits(session_files_subj{isesh});
            if any(isnan([rich_waits; poor_waits])) && isesh==1 
                try
                [rich_waits, poor_waits] = rich_poor_waits(session_files_subj{isesh+1});
                catch
                end
            end
            if any(isnan([rich_waits; poor_waits])) && isesh==length(session_files_subj) 
                try
                [rich_waits, poor_waits] = rich_poor_waits(session_files_subj{isesh-1});
                catch
                end
            end

            %load session zdiff
            total_learn_NL(isubj_NL, sesh_ct) = zdiff(rich_waits, poor_waits);

        end
    end

    % delta zdiff
    total_learn_NL = total_learn_NL(:,2)-total_learn_NL(:,1);
    total_newlearn(inl,:) = [nanmean(total_learn_NL,1) nanstd(total_learn_NL,[],1)./sqrt(sum(~isnan(total_learn_NL),1))];
end



% plot
%
figure; hold on
for inl = 1:7
    % training
    errorbar(1:length(total_training(:,1,inl)), total_training(:,1,inl), total_training(:,2,inl), '-', 'markersize', 20, 'color', colors(inl,:))
end
for inl = 1:7
    % new learning
    errorbar(num_training_stages+1+inl/2, total_newlearn(inl,1), total_newlearn(inl,2), '.', 'markersize', 20, 'color', colors(inl,:))
end
ylabel('Total learning (z)')
xlabel('Training Stage')
set(gca,'TickLength',[0, 0]); box off;
xlim([0 16])
ylim([-.5 5])
xticks([1:num_training_stages])
plot(xlim, [1 1].*0, 'k--')
legend({'NL01','NL02','NL03','NL04','NL05','NL06','NL07'}, 'location', 'northeastoutside')

end

function [rich_waits, poor_waits] = rich_poor_waits(session_file_string)
            
% load session
load(session_file_string, 'trl_mtx', 'medass_cell');

% wait times
[wait_durations, aw_freq] = wait_times_prep(trl_mtx,2);
[prob_dist, pd_freq] = rwd_prob_by_freq(medass_cell);

% rich and poor waits
min_samps = 2;
rich_waits = wait_durations(ismember(aw_freq,pd_freq(prob_dist>0.5)));       
if length(rich_waits) < min_samps
    rich_waits = nan; 
end
poor_waits = wait_durations(ismember(aw_freq,pd_freq(prob_dist<0.5)));
if length(poor_waits) < min_samps
    poor_waits = nan; 
end

end

function B=nancumsum(A,dim,nmode)
%
if nargin < 3
    nmode = 1;
end
if ~ismember(nmode,1:4)
    error('NANCUMSUM: unacceptable value for nmode parameter.');
end
if nargin < 2 || isempty(dim)
    if ~isscalar(A)
        dim = find(size(A)>1);
        dim = dim(1);
    else
        % For scalar inputs (no nonsingleton dimension)
        dim = 1;
    end
end
% Calculate cumulative sum, depending on selection of nmode
switch nmode
    case 1
        % TREAT NaNs as 0's
        B = A;
        B(B~=B) = 0;
        B = cumsum(B, dim);
    case 2
        % DO NOT INCREMENT, BUT USE NaNs AS PLACEHOLDERS.
        B = nancumsum(A, dim, 1);
        B(A~=A) = NaN;
     case 3
        % RESET sum on NaNs, replacing NaNs with zeros.
        naninds = find(A~=A);
        for ii = 1:numel(naninds)
            B = nancumsum(A, dim, 1);
            A(naninds(ii)) = -B(naninds(ii));
        end
        B = cumsum(A,dim);
    otherwise %case 4
        % RESET sum on NaNs, maintaining NaNs as position holders.
        naninds = find(A~=A);
        B = nancumsum(A,dim,3);
        B(naninds)= NaN;
end

end

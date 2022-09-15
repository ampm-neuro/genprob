function ALL_opto_discrim(file_keywords, first_last)
% combine cells from a single stage from all animals

% get all files
fp = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_optoExp_hpc\';%train_mevar_optoExp_acc\';
session_files_mevar = get_file_paths_targeted(fp, file_keywords);
fp = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_hivar_optoExp_hpc\';%train_mevar_optoExp_acc\';
session_files_hivar = get_file_paths_targeted(fp, file_keywords);

session_files = [session_files_mevar; session_files_hivar];

% find unique subjects
subjs = [];
for isf = 1:size(session_files,1)
    subjs = [subjs; {session_files{isf}(strfind(session_files{isf}, '\\')+2 : strfind(session_files{isf}, '\\')+9)}];
end
subjs = unique(subjs);

% find one file per subject
session_files_hold = [];
for isubj = 1:length(subjs)
    subj_sf = session_files(contains(session_files, subjs{isubj}));
    if strcmp(first_last, 'first')
        session_files_hold = [session_files_hold; subj_sf(1)];
    elseif strcmp(first_last, 'last')
        if size(subj_sf,1) == 1
            continue
        end
        session_files_hold = [session_files_hold; subj_sf(end)];
    elseif strcmp(first_last, 'all')
        session_files_hold = [session_files_hold; subj_sf(:)];
    else
        error('do better with first_last input')
    end
end
session_files = session_files_hold


% iterate through session files
all_means = nan(size(session_files,1), 4); % poor nolas, rich nolas, poor las, rich las
all_freqs = nan(size(session_files,1), 4); % poor, rich
for isf = 1:size(session_files,1)

    % load session file
    load(session_files{isf}, 'medass_cell', 'trl_mtx', 'unq_frq')
        
    % idx
    probe_idx = trl_mtx(:,3)==0;
    optoON_idx = trl_mtx(:,13)==1;
    rich_idx = rich_trl_idx(trl_mtx); 
    
    % load
    all_means(isf, 1) = mean(trl_mtx(probe_idx & optoON_idx & ~rich_idx, 12));
    all_means(isf, 2) = mean(trl_mtx(probe_idx & optoON_idx & rich_idx, 12));
    all_means(isf, 3) = mean(trl_mtx(probe_idx & ~optoON_idx & ~rich_idx, 12));
    all_means(isf, 4) = mean(trl_mtx(probe_idx & ~optoON_idx & rich_idx, 12));
    all_freqs(isf, [1 3]) = mode(floor(trl_mtx(~rich_idx, 2)));
    all_freqs(isf, [2 4]) = mode(floor(trl_mtx(rich_idx, 2)));
    
    % clear lingering variables
    clearvars('medass_cell', 'trl_mtx', 'unq_frq')
end

% remove incomplete sessions
all_means(isnan(sum(all_means,2))) = nan;


% plot colors 
current_colors = [0 0 1; .5 .5 .5]; %distinguishable_colors(2);

% plot dots
jxpos = nan(size(all_freqs)); pct_jit = 7;
figure;
hold on
legend_trick(current_colors, '-')
for icol = 1:2
     jxpos(:,icol) = all_freqs(:, icol);
     jxpos(:,icol) = ((randi([100-pct_jit,100+pct_jit], [size(jxpos,1),1]).*jxpos(:,icol))./100);
     plot(jxpos(:,icol), all_means(:, icol), 'o', 'color', current_colors(1,:))
     
     jxpos(:,icol+2) = all_freqs(:, icol+2);
     jxpos(:,icol+2) = ((randi([100-pct_jit,100+pct_jit], [size(jxpos,1),1]).*jxpos(:,icol+2))./100);
     plot(jxpos(:,icol+2), all_means(:, icol+2), 'o', 'color', current_colors(2,:))
end


% plot lines connecting dots
for isubj = 1:size(all_means,1)
    plot(jxpos(isubj,[1 2]), all_means(isubj, [1 2]), '-', 'color', current_colors(1,:), 'linewidth', 0.5)
    plot(jxpos(isubj,[3 4]), all_means(isubj, [3 4]), '-', 'color', current_colors(2,:), 'linewidth', 0.5)
end

% plot errorbars
errorbar(mode(all_freqs(:,[1 2])), nanmean(all_means(:,[1 2])), nanstd(all_means(:,[1 2]))./sqrt(size(all_means,1)),...
    '-', 'color', current_colors(1,:), 'linewidth', 2.5);
errorbar(mode(all_freqs(:,[3 4])), nanmean(all_means(:,[3 4])), nanstd(all_means(:,[3 4])./sqrt(size(all_means,1))),...
    '-', 'color', current_colors(2,:), 'linewidth', 2.5);

% aesthetics
set(gca,'TickLength',[0, 0]);
xlabel('Tone Frequency (Hz)')
ylabel('Wait Durations (s)')
ylim_hold = ylim; ylim([0 ylim_hold(2)]);
set(gca, 'XScale', 'log')
title('Waits')
xlim([4500 42000])
xticks([5000 8500 14000 23000 35000])
legend({'LaserON', 'LaserOFF'})



% cellCorrs_v_probeCorrs
% compares average correlation between shared cells with correlation
% between probe response profiles


% subject ids
subject_ids = [{'651049m1'} {'658648m2'} {'691359m1'} {'690330m1'} {'690330m2'} {'691359m2'} {'696944m3'} {'699437m4'} {'699438m3'}...
    {'683472m2'} {'683472m3'} {'687034m4'} {'699437'}]
hivar_mevar = [{'mevar'} {'mevar'} {'mevar'} {'mevar'} {'mevar'} {'mevar'} {'mevar'} {'mevar'} {'mevar'}...
    {'hivar'} {'hivar'} {'hivar'} {'hivar'}];


%[{'651049m1'} {'658648m2'} {'683472m2'} {'683472m3'}];
%hivar_mevar = [{'mevar'} {'mevar'} {'hivar'} {'hivar'}];

%% compute correlation between all cells all sessions
%
 [all_cell_corrs_mevar, all_merge_mtx_mevar, subj_corr_cell_mevar] = cell_turnover_timewarp_trials_multisubj('mevar', subject_ids([1:5]));
 [all_cell_corrs_hivar, all_merge_mtx_hivar, subj_corr_cell_hivar] = cell_turnover_timewarp_trials_multisubj('hivar', subject_ids([6:8]));
 save('cellCorrs_v_probeCorrs_prep_2.mat')
 %}
%load('cellCorrs_v_probeCorrs_prep_2.mat')

% probe days only
probe_days = [1:3:19 20];
cell_corrs = [];
for imevar = 1:sum(strcmp(hivar_mevar, 'mevar'))
    cell_corrs = [cell_corrs {subj_corr_cell_mevar{imevar}(probe_days(probe_days<=length(subj_corr_cell_mevar{imevar})), probe_days(probe_days<=length(subj_corr_cell_mevar{imevar})))}];
end
for ihivar = 1:sum(strcmp(hivar_mevar, 'hivar'))
    cell_corrs = [cell_corrs {subj_corr_cell_hivar{ihivar}(probe_days(probe_days<=length(subj_corr_cell_hivar{ihivar})), probe_days(probe_days<=length(subj_corr_cell_hivar{ihivar})))}];
end
%cell_corrs = [{subj_corr_cell_mevar{1}(probe_days(probe_days<=length(subj_corr_cell_mevar{1})), probe_days(probe_days<=length(subj_corr_cell_mevar{1})))}...
%                {subj_corr_cell_mevar{2}(probe_days(probe_days<=length(subj_corr_cell_mevar{2})), probe_days(probe_days<=length(subj_corr_cell_mevar{2})))}...
%                {subj_corr_cell_hivar{1}(probe_days(probe_days<=length(subj_corr_cell_hivar{1})), probe_days(probe_days<=length(subj_corr_cell_hivar{1})))}...
%                {subj_corr_cell_hivar{2}(probe_days(probe_days<=length(subj_corr_cell_hivar{2})), probe_days(probe_days<=length(subj_corr_cell_hivar{2})))}];

            
% minimum number of cells
min_cell = 5;

% means only
cell_corr_means = nan(8,8,length(cell_corrs));
for isubj = 1:length(cell_corrs)
    for isesh1 = 1:size(cell_corrs{isubj},1)
        for isesh2 = 1:size(cell_corrs{isubj},2)
            
            % no redundancy
            if isesh2>=isesh1
                continue
            end
            if length(cell_corrs{isubj}{isesh1, isesh2})>=min_cell
                cell_corr_means(isesh1,isesh2,isubj) = nanmean(cell_corrs{isubj}{isesh1, isesh2});
            else
                cell_corr_means(isesh1,isesh2,isubj) = nan;
            end
        end
    end
end



%% get probe wait times
[~, subj_ids_mevar, all_waits_mevar] = probe_wait_times('train_mevar_imaging_hpc_fin', 1:8, 1:41);
[~, subj_ids_hivar, all_waits_hivar] = probe_wait_times('train_hivar_imaging_hpc_fin', 1:8, 1:41);

% subj id order should be same as above
subj_ids_mevar = subj_ids_mevar{1};
subj_ids_hivar = subj_ids_hivar{1};
subj_id_probes = [subj_ids_mevar;subj_ids_hivar];

% merge probe wait times
all_waits_probes = cell(1,8);
for iprobe = 1:8
    all_waits_probes{iprobe} = [all_waits_mevar{iprobe};all_waits_hivar{iprobe}];
end

% compute probe corms
probe_corrs = nan(8,8,length(subject_ids));
for isubj = 1:length(subject_ids)
    for iprobe1 = 1:length(all_waits_probes)
        for iprobe2 = 1:length(all_waits_probes)
            
            % no redundancy
            if iprobe2>=iprobe1
                continue
            end
            
            % waits
            probe1 = all_waits_probes{iprobe1}(isubj,:)';
            probe2 = all_waits_probes{iprobe2}(isubj,:)';
            
            % load corrs
            probe_corrs(iprobe1,iprobe2,isubj) = corr(probe1, probe2);
        end
    end
end 

%% plot

% figure
figure; hold on
% plot circles
for isubj = 1:length(subject_ids)
    
    try
    
    % inputs
    probe_corrs_local = probe_corrs(:,:,isubj);
    cell_corr_means_local = cell_corr_means(:,:,isubj);
    
    cell_corr_means_local = cell_corr_means_local(:);
    probe_corrs_local = probe_corrs_local(:);
    
    cell_corr_means_local = atanh(cell_corr_means_local);
    probe_corrs_local = atanh(probe_corrs_local);
    
    % color
    if strcmp(hivar_mevar{isubj},'mevar')
        plot(cell_corr_means_local(:), probe_corrs_local(:), 'go')
    elseif strcmp(hivar_mevar{isubj},'hivar')
        plot(cell_corr_means_local(:), probe_corrs_local(:), 'bo')
    else
        error('mevar hivar input issue')
    end
    
    catch
    end
end


% plot fit line

cell_corr_means = cell_corr_means(:);
probe_corrs = probe_corrs(:);

cell_corr_means = atanh(cell_corr_means);
probe_corrs = atanh(probe_corrs);

%probe_corrs(probe_corrs<-.5) = nan;
%probe_corrs(probe_corrs>.6) = nan;
[r,p] = fit_line(cell_corr_means, probe_corrs,0);

title(['r=' num2str(r) '; p=' num2str(round(p*10000)/10000)])
xlabel('Population firing similarity (r)')
ylabel('Memory similarity (r)') 

tic_vect = [-.9 -.8 -.6 -.3 0 .3 .6 .8 .9]; 
xlim(atanh([min(tic_vect) max(tic_vect)])); xticks(atanh(tic_vect)); xticklabels(tic_vect)
ylim(atanh([min(tic_vect) max(tic_vect)])); yticks(atanh(tic_vect)); yticklabels(tic_vect)

xlim(atanh([-.9 .9])); xticks(atanh(tic_vect)); xticklabels(tic_vect)
ylim(atanh([-.9 .9])); yticks(atanh(tic_vect)); yticklabels(tic_vect)

hold on; plot([0 0], ylim, 'k--')
hold on; plot(xlim, [0 0], 'k--')

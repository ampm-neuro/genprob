% PAIRWISE overlap population vs probe similarity

%
%% subject ids
subject_ids = [{'651049m1'} {'658648m2'} {'691359m1'} {'690330m1'} {'690330m2'} {'691359m2'} {'696944m3'} {'699437m4'} {'699438m3'}...
    {'683472m2'} {'683472m3'} {'687034m4'} {'699437m2'} {'699438m1'}];
hivar_mevar = [{'mevar'} {'mevar'} {'mevar'} {'mevar'} {'mevar'} {'mevar'} {'mevar'} {'mevar'} {'mevar'}...
    {'hivar'} {'hivar'} {'hivar'} {'hivar'} {'hivar'}];

%
%% population overlap probes only
probe_idx = [1:3:19 20];

% mevar
overlap_mtx_mevar = overlapping_population_average(subject_ids(strcmp(hivar_mevar, 'mevar')));
overlap_mtx_probes_mevar = overlap_mtx_mevar(probe_idx, probe_idx, :); 

% hivar
overlap_mtx_hivar = overlapping_population_average(subject_ids(strcmp(hivar_mevar, 'hivar')));
overlap_mtx_probes_hivar = overlap_mtx_hivar(probe_idx, probe_idx, :); 

% merge
overlap_mtx_probes = cat(3, overlap_mtx_probes_mevar, overlap_mtx_probes_hivar);
%


%% full pairwise loop
for isubj = 1:size(overlap_mtx_probes,3)
    for iprobe1 = 1:size(overlap_mtx_probes,1)
        for iprobe2 = 1:size(overlap_mtx_probes,2)
            if iprobe2>=iprobe1
                
                % population overlap
                overlap_mtx_probes(iprobe1, iprobe2, isubj) = nan;
                
            end
        end
    end
end


%% get probe wait time similarity
[~, subj_ids_mevar, all_waits_mevar] = probe_wait_times('train_mevar_imaging_hpc_fin', 1:8, 1:41);
[~, subj_ids_hivar, all_waits_hivar] = probe_wait_times('train_hivar_imaging_hpc_fin', 1:8, 1:41);

size(all_waits_mevar)
size(all_waits_hivar)

% subj id order should be same as above
subj_ids_mevar = subj_ids_mevar{1};
subj_ids_hivar = subj_ids_hivar{1};
subj_id_probes = [subj_ids_mevar;subj_ids_hivar]

% merge probe wait times
all_waits_probes = cell(1,8);
for iprobe = 1:8
    all_waits_probes{iprobe} = [all_waits_mevar{iprobe}; all_waits_hivar{iprobe}];
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
            
            iprobe1
            iprobe2
            isubj
            
            % waits
            probe1 = all_waits_probes{iprobe1}(isubj,:)';
            probe2 = all_waits_probes{iprobe2}(isubj,:)';
            
            % load corrs
            probe_corrs(iprobe1,iprobe2,isubj) = corr(probe1, probe2);
        end
    end
end 



%}

%% plot

% figure
figure; hold on

% plot circles
for isubj = 1:length(subject_ids)

        % inputs
        probe_sim_local = probe_corrs(:,:,isubj);
        pop_overlap_local = overlap_mtx_probes(:,:,isubj);

        pop_overlap_local = pop_overlap_local(:);
        probe_sim_local = probe_sim_local(:);

        % color
                
        if strcmp(hivar_mevar{isubj},'mevar')
            plot(pop_overlap_local(:), probe_sim_local(:), 'go')
        elseif strcmp(hivar_mevar{isubj},'hivar')
            plot(pop_overlap_local(:), probe_sim_local(:), 'bo')
        else
            error('mevar hivar input issue')
        end

end


pop_overlap = overlap_mtx_probes(:);
probe_sim = probe_corrs(:);

% atanh r values
probe_sim = atanh(probe_sim);

[r,p] = fit_line(pop_overlap, probe_sim,0);

title(['r=' num2str(r) '; p=' num2str(round(p*10000)/10000)])
xlabel('Proportion pop overlap')
ylabel('Memory similarity (r)') 

tic_vect = [-.9 -.8 -.6 -.3 0 .3 .6 .8 .9]; 
xlim(atanh([min(tic_vect) max(tic_vect)])); xticks(atanh(tic_vect)); xticklabels(tic_vect)
xlim(atanh([-.9 .9])); xticks(atanh(tic_vect)); xticklabels(tic_vect)
ylim(atanh([min(tic_vect) max(tic_vect)])); yticks(atanh(tic_vect)); yticklabels(tic_vect)
ylim(atanh([-.9 .9])); yticks(atanh(tic_vect)); yticklabels(tic_vect)

hold on; plot([0 0], ylim, 'k--')
hold on; plot(xlim, [0 0], 'k--')


%% mevar hivar individual plots

% mevar
overlap_vect = colon_op(overlap_mtx_probes(:,:,strcmp(hivar_mevar, 'mevar')));
probe_vect = colon_op(probe_corrs(:,:,strcmp(hivar_mevar, 'mevar')));

%overlap_vect(probe_vect<-.4)=[];
%probe_vect(probe_vect<-.4)=[];

figure; [r,p] = fit_line(overlap_vect , atanh(probe_vect)); 
title(['mevar; r=' num2str(r) '; p=' num2str(round(p*10000)/10000)])
xlabel('Proportion pop overlap')
ylabel('Memory similarity (r)') 
ylim(atanh([min(tic_vect) max(tic_vect)])); yticks(atanh(tic_vect)); yticklabels(tic_vect)
ylim(atanh([-.6 .6])); yticks(atanh(tic_vect)); yticklabels(tic_vect)

hold on; plot([0 0], ylim, 'k--')
hold on; plot(xlim, [0 0], 'k--')


% hivar
overlap_vect = colon_op(overlap_mtx_probes(:,:,strcmp(hivar_mevar, 'hivar')));
probe_vect = colon_op(probe_corrs(:,:,strcmp(hivar_mevar, 'hivar')));

figure; [r,p] = fit_line(overlap_vect , atanh(probe_vect)); 
title(['hivar; r=' num2str(r) '; p=' num2str(round(p*10000)/10000)])
xlabel('Proportion pop overlap')
ylabel('Memory similarity (r)') 
ylim(atanh([min(tic_vect) max(tic_vect)])); yticks(atanh(tic_vect)); yticklabels(tic_vect)
ylim(atanh([-.6 .6])); yticks(atanh(tic_vect)); yticklabels(tic_vect)

hold on; plot([0 0], ylim, 'k--')
hold on; plot(xlim, [0 0], 'k--')




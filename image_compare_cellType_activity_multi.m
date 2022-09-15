function [all_dm, all_rm, all_rc, eb_in] = image_compare_cellType_activity_multi(mevar_or_hivar, subject_ids)
% plot reliability as a function of cell type, reactivated or new


% ('mevar', [{'651049m1'} {'658648m2'}]);
%('hivar', [{'683472m2'} {'683472m3'}]);

% chron
tses = [0.2 1.1 1.0 2.0 8.5 2.0];
session_chron_reorder_idx = [1 9 15 2 10 16 3 11 17 4 12 18 5 13 19 6 14 20 7 8];

%% session cell
all_session_cell = cell(1,length(subject_ids));
all_crm = cell(1,length(subject_ids));
for isubj = 1:length(subject_ids)

    
    % load cell registration matrix
    load(['cell_reg_' subject_ids{isubj} '.mat'])
    cell_regist_mtx = cell_registered_struct.cell_to_index_map;
    
    
    % chron
    if size(cell_regist_mtx,2)==19
        scri = session_chron_reorder_idx(1:end-1); scri(scri>=8) = scri(scri>=8)-1;
    else
        scri = session_chron_reorder_idx;
    end
    
    all_crm{isubj} = cell_regist_mtx(:, scri);
    
    % get session paths
    last_prob_sessions = [];
    for iprob = 1:6
        prob_sessions = get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_' mevar_or_hivar '_imaging_hpc\' subject_ids{isubj} ], ['var0' num2str(iprob)], 'LED');
        last_prob_sessions = [last_prob_sessions; prob_sessions(end)];
    end
    session_cell = [...
        get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_' mevar_or_hivar '_imaging_hpc\' subject_ids{isubj} ], 'preprobe', 'LED');...
        get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_' mevar_or_hivar '_imaging_hpc\' subject_ids{isubj} ] , 'postprobe', 'LED');...
        get_file_paths_targeted(['C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_' mevar_or_hivar '_imaging_hpc\' subject_ids{isubj} ], 'var0', '-01', 'LED');...
        last_prob_sessions];
    all_session_cell{isubj} = session_cell(scri);
    
end






%% compute
% for each session
all_dm = [];
all_rm = [];
all_rc = [];
for isubj = 1:length(subject_ids)
    [dispersion_mtx, reliability_mtx, reactivation_ct_mtx] = image_compare_cellType_activity_reactct(all_session_cell{isubj}, all_crm{isubj}, tses);
    
    % create dummies to catch subjects within less than 20 sessions
    dispersion_mtx_hold = nan(size(dispersion_mtx,1), 20);
    reliability_mtx_hold = nan(size(reliability_mtx,1), 20);
    reactivation_ct_mtx_hold = nan(size(reactivation_ct_mtx,1), 20);
    
    % load dummies
    dispersion_mtx_hold(:, 1:size(dispersion_mtx,2)) = dispersion_mtx;
    reliability_mtx_hold(:, 1:size(reliability_mtx,2)) = reliability_mtx;
    reactivation_ct_mtx_hold(:, 1:size(reactivation_ct_mtx,2)) = reactivation_ct_mtx;
    
    all_dm = [all_dm; dispersion_mtx_hold];
    all_rm = [all_rm; reliability_mtx_hold];
    all_rc = [all_rc; reactivation_ct_mtx_hold];
end


%% plot

eb_in = cell(1,max(all_rc(:)));
try
for irc = 1:max(all_rc(:))
    eb_in{irc} = all_rm(all_rc==irc);
end
figure; errorbar_plot(eb_in(1:10))

% cor mtx
mtx = nan(10,20); 
for isesh = 1:20
    rm = all_rc(:,isesh);  
    % max 10 reacts
    for ict = 1:10
        if ict<10
            meanr = mean(all_rm(rm==ict,isesh)); 
        else
            meanr = mean(all_rm(rm>=ict,isesh)); 
        end
        mtx(ict, isesh) = meanr; 
    end
end
figure; imagesc(mtx)
catch
end
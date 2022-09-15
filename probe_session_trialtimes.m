%probe_session_trialtimes
% get median duration of each trial event

% all session files
sesh_files_mevar = get_file_paths_targeted('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_mevar_imaging_hpc\', '.mat');
sesh_files_mevar = sesh_files_mevar(contains(sesh_files_mevar, 'probe_0'));
sesh_files_hivar = get_file_paths_targeted('C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\train_hivar_imaging_hpc\', '.mat');
sesh_files_hivar = sesh_files_hivar(contains(sesh_files_mevar, 'probe_0'));

sesh_files = [sesh_files_mevar; sesh_files_hivar];

[all_trl_mtx, session_number, problem_number] = ALL_trl_mtx(sesh_files);

median_nosepoke = median(all_trl_mtx(:,9)-all_trl_mtx(:,6));
median_nosepoke_to_headentry = median(all_trl_mtx(:,9)-all_trl_mtx(:,6));
median_delay = median(all_trl_mtx(:,12));

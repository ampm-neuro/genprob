function save_footprints_cellregprep(mouse, folder)
% save footprints in each mouse folder

folder = 'train_hivar_imaging_hpc';
mouse = '699437m2'; %'658648m2';

cell_reg_temp_path = ['C:\Users\ampm1\Documents\MATLAB\generalization\cellreg_data\' mouse '\spatial_footprints'];

% sessions
%
for isesh = 7:12
    
    isesh
    
    prob_num = isesh-6
    
    if contains(folder, 'mevar')
    
        %{
        if isesh < 10
            folder_path = ['E:\two_tone\' folder '\' mouse '\LED_gen0' num2str(isesh) '_mevar0' num2str(prob_num)];
        else
            folder_path = ['E:\two_tone\' folder '\' mouse '\LED_gen' num2str(isesh) '_mevar0' num2str(prob_num)];
        end
        %}
        
        if isesh < 10
            if exist(['E:\two_tone\' folder '\' mouse '\LED_gen0' num2str(isesh) '_mevar0' num2str(prob_num)], 'file')
            	folder_path = ['E:\two_tone\' folder '\' mouse '\LED_gen0' num2str(isesh) '_mevar0' num2str(prob_num)];
            else
                folder_path = ['E:\two_tone\' folder '\' mouse '\LED_gen0' num2str(isesh) '_mevar0' num2str(prob_num) '_rev'];
            end
        else
            if exist(['E:\two_tone\' folder '\' mouse '\LED_gen' num2str(isesh) '_mevar0' num2str(prob_num)], 'file')
                folder_path = ['E:\two_tone\' folder '\' mouse '\LED_gen' num2str(isesh) '_mevar0' num2str(prob_num)];
            else
                folder_path = ['E:\two_tone\' folder '\' mouse '\LED_gen' num2str(isesh) '_mevar0' num2str(prob_num) '_rev'];
            end
        end
    
    elseif contains(folder, 'hivar')
        
        if isesh < 10
            if exist(['E:\two_tone\' folder '\' mouse '\LED_gen0' num2str(isesh) '_hivar0' num2str(prob_num)], 'file')
            	folder_path = ['E:\two_tone\' folder '\' mouse '\LED_gen0' num2str(isesh) '_hivar0' num2str(prob_num)];
            else
                folder_path = ['E:\two_tone\' folder '\' mouse '\LED_gen0' num2str(isesh) '_hivar0' num2str(prob_num) '_rev'];
            end
        else
            if exist(['E:\two_tone\' folder '\' mouse '\LED_gen' num2str(isesh) '_hivar0' num2str(prob_num)], 'file')
                folder_path = ['E:\two_tone\' folder '\' mouse '\LED_gen' num2str(isesh) '_hivar0' num2str(prob_num)];
            else
                folder_path = ['E:\two_tone\' folder '\' mouse '\LED_gen' num2str(isesh) '_hivar0' num2str(prob_num) '_rev'];
            end
        end
        
    else
        error
    end
        
    session_folders = get_folder_paths_all(folder_path);
    session_paths = get_file_paths_targeted(folder_path, '.Subject ');
    
    for first_last = 1:2
        
        if first_last ==1
            session_folder = session_folders{1}; session_folder = session_folder(1:strfind(session_folder, '\01')-1);

            try
            save_footprints(session_folders{1});
            
            
            
            %[session_folders{1} '\spatial_footprints.mat'], [cell_reg_temp_path '\spatial_footprints_problem_first_0' num2str(isesh-6) '.mat']
            copyfile([session_folders{1} '\spatial_footprints.mat'], [cell_reg_temp_path '\spatial_footprints_problem_first_0' num2str(isesh-6) '.mat'])
            
            save_footprints([session_folder '\01']);
            
            
            catch
            %    session_folders{1}
            end
            
            
            %{
            if isesh < 10
                save_footprints([session_folder '\01']);
            else
                save_footprints([session_folder '\01']);
            end
            %}
        else
            save_footprints(session_folders{end});
            [session_folders{end} '\spatial_footprints.mat'], [cell_reg_temp_path '\spatial_footprints_problem_last_0' num2str(isesh-6) '.mat']
            
            [cell_reg_temp_path '\spatial_footprints_problem_last_0' num2str(isesh-6) '.mat']
            
            copyfile([session_folders{end} '\spatial_footprints.mat'], [cell_reg_temp_path '\spatial_footprints_problem_last_0' num2str(isesh-6) '.mat'])
            %
            if size(session_paths,1) < 10
                save_footprints([session_folder '\0' num2str(size(session_paths,1))]);
            else
                save_footprints([session_folder '\' num2str(size(session_paths,1))]);
            end
            %
        end
        
    end
end
%}
% probes
for iprobe = 1:8
    iprobe
    if exist(['E:\two_tone\' folder '\' mouse '\LED_gen14_probe\0' num2str(iprobe)], 'file')
        save_footprints(['E:\two_tone\' folder '\' mouse '\LED_gen14_probe\0' num2str(iprobe)])
        copyfile(['E:\two_tone\' folder '\' mouse '\LED_gen14_probe\0' num2str(iprobe) '\spatial_footprints.mat'], [cell_reg_temp_path '\spatial_footprints_probe0' num2str(iprobe) '.mat'])
    else
        save_footprints(['E:\two_tone\' folder '\' mouse '\LED_gen14_probe_rev\0' num2str(iprobe)])
    	copyfile(['E:\two_tone\' folder '\' mouse '\LED_gen14_probe_rev\0' num2str(iprobe) '\spatial_footprints.mat'], [cell_reg_temp_path '\spatial_footprints_probe0' num2str(iprobe) '.mat'])
    end
end




function load_medass_all(fpath)
% runs load_medass on all files in a folder
% saves in individual subject folders at 'destp'
%   place to-be-filed files into correct experimental folder, then this
%   program will place into correct subject, stage, and day folders

%destination path (contains all experimental folders)
destp = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data';

%all files

%fpath
origin_paths_all = get_file_paths_all(fpath);
origin_paths_all = origin_paths_all(~contains(origin_paths_all, 'old'));
origin_paths_all = origin_paths_all(~contains(origin_paths_all, '.mat'));

probe_count = [];
for ifile = 1:length(origin_paths_all)
    
    %path folders
    path_segments = regexp(origin_paths_all{ifile}, filesep, 'split');

    %relevant folder information
    task = path_segments{end-5};
    sub_folder = path_segments{end-4};
    subject_id = path_segments{end-3};
    stage = path_segments{end-2};
        % remove rev from stage name
        if contains(stage, '_rev')
            stage(strfind(stage, '_rev') : strfind(stage, '_rev') +3) = [];
        end
    sesh = path_segments{end-1};
    
    %ensure required folders exist
    while true
        if strcmp(task, sub_folder)
            if ~exist(destp, 'dir')
                mkdir(destp)
                continue
            elseif ~exist([destp '\' task] , 'dir')
                mkdir([destp '\' task])
                continue
            elseif ~exist([destp '\' task '\' subject_id] , 'dir')
                mkdir([destp '\' task '\' subject_id])
                continue
            end
            path_hold = [destp '\' task '\' subject_id];
            break
        else
            if ~exist(destp, 'dir')
                mkdir(destp)
                continue
            elseif ~exist([destp '\' task] , 'dir')
                mkdir([destp '\' task])
                continue
            elseif ~exist([destp '\' task '\' sub_folder] , 'dir')
                mkdir([destp '\' task '\' sub_folder])
                continue
            elseif ~exist([destp '\' task '\' sub_folder '\' subject_id] , 'dir')
                mkdir([destp '\' task '\' sub_folder '\' subject_id])
                continue
            end
            path_hold = [destp '\' task '\' sub_folder '\' subject_id];
            break
        end
    end

    %saved filed name
    save_file_name = [stage '-' sesh];
    
    
    
    %% probe
    if contains(save_file_name, 'probe') && ~contains(save_file_name, 'quiet')
        
        % info from current file
        [~, current_train_date, subject_id] = medass_mpc_id(origin_paths_all{ifile});
        
        % all origin paths containing subject id
        origin_paths_subj = origin_paths_all(contains(origin_paths_all, subject_id));
        
        % training days only       
        dest_paths_subj = get_file_paths_all(path_hold); % destination folder
        dest_paths_subj_train = dest_paths_subj(contains(dest_paths_subj, {'novar0', 'lovar0', 'mevar0', 'hivar0', 'exvar0', 'ctl0'}));
        origin_paths_subj_train = origin_paths_subj(contains(origin_paths_subj, {'novar0', 'lovar0', 'mevar0', 'hivar0', 'exvar0', 'ctl0'}));
        
        % iterate through destination ids finding corresponding origin
        % files
        comb_dest_origin_idx = false(size(origin_paths_subj_train,1),1);
        for idids = 1:size(dest_paths_subj_train,1)
            
            % training id
            dest_ids_1 = [dest_paths_subj_train{idids}(strfind(dest_paths_subj_train{idids},'novar0'):strfind(dest_paths_subj_train{idids},'novar0')+6) ...
                dest_paths_subj_train{idids}(strfind(dest_paths_subj_train{idids},'lovar0'):strfind(dest_paths_subj_train{idids},'lovar0')+6) ...
                dest_paths_subj_train{idids}(strfind(dest_paths_subj_train{idids},'mevar0'):strfind(dest_paths_subj_train{idids},'mevar0')+6) ...
                dest_paths_subj_train{idids}(strfind(dest_paths_subj_train{idids},'hivar0'):strfind(dest_paths_subj_train{idids},'hivar0')+6) ...
                dest_paths_subj_train{idids}(strfind(dest_paths_subj_train{idids},'exvar0'):strfind(dest_paths_subj_train{idids},'exvar0')+6) ...
                dest_paths_subj_train{idids}(strfind(dest_paths_subj_train{idids},'ctl0'):strfind(dest_paths_subj_train{idids},'ctl0')+4)];
            % session id
            dest_ids_2 = dest_paths_subj_train{idids}(strfind(dest_paths_subj_train{idids},'-')+1:strfind(dest_paths_subj_train{idids},'-')+2);

            % find corresponding origin files            
            if contains(dest_paths_subj_train{idids}, '_opto-')
                comb_dest_origin_idx = comb_dest_origin_idx | contains(origin_paths_subj_train, [dest_ids_1 '_opto\' dest_ids_2])...
                    | contains(origin_paths_subj_train, [dest_ids_1 '_rev_opto\' dest_ids_2]);
            else
                comb_dest_origin_idx = comb_dest_origin_idx | contains(origin_paths_subj_train, [dest_ids_1 '\' dest_ids_2])...
                    | contains(origin_paths_subj_train, [dest_ids_1 '_4t\' dest_ids_2])...
                    | contains(origin_paths_subj_train, [dest_ids_1 '_2t\' dest_ids_2])...
                    | contains(origin_paths_subj_train, [dest_ids_1 '_rev\' dest_ids_2]);
            end
            
        end
        origin_paths_subj_train_indest = origin_paths_subj_train(comb_dest_origin_idx);
        
        if size(dest_paths_subj_train,1) ~= size(origin_paths_subj_train_indest,1) 
            origin_paths_subj_train_indest
            dest_paths_subj_train
            %
            if size(dest_paths_subj_train,1) > size(origin_paths_subj_train_indest,1) 
                uncommon_paths = setdiff(origin_paths_subj_train_indest, dest_paths_subj_train)
            else
                uncommon_paths = setdiff(dest_paths_subj_train, origin_paths_subj_train_indest)
            end
            %}
            error('origin and destination file list mismatch')
        end
        
        % training data files dated before the current probe file
        path_dates = cell(size(origin_paths_subj_train_indest,1),1);
        for ipd = 1:size(origin_paths_subj_train_indest,1)
            [~,path_dates{ipd}] = medass_mpc_id(origin_paths_subj_train_indest{ipd});
            fclose('all');
        end
                
        if ~isempty(path_dates) && ~isempty(dest_paths_subj_train(datenum(path_dates)<datenum(current_train_date)))
            
            % only dates before current session
            dest_paths_subj_train = dest_paths_subj_train(datenum(path_dates)<datenum(current_train_date));
            origin_paths_subj_train_indest = origin_paths_subj_train_indest(datenum(path_dates)<datenum(current_train_date));
            
            % search origin folder using info from destination folder
            last_training_path = origin_paths_subj_train_indest{end};
            [~, last_train_date] = medass_mpc_id(last_training_path);

            if contains(last_training_path, 'novar0')
                last_train_num = str2double(last_training_path(strfind(last_training_path, 'novar0')+6));
            elseif  contains(last_training_path, 'lovar0')
                last_train_num = str2double(last_training_path(strfind(last_training_path, 'lovar0')+6));
            elseif  contains(last_training_path, 'mevar0')
                last_train_num = str2double(last_training_path(strfind(last_training_path, 'mevar0')+6));
            elseif  contains(last_training_path, 'hivar0')
                last_train_num = str2double(last_training_path(strfind(last_training_path, 'hivar0')+6));
            elseif  contains(last_training_path, 'exvar0')
                last_train_num = str2double(last_training_path(strfind(last_training_path, 'exvar0')+6));
            elseif  contains(last_training_path, 'ctl0')
                last_train_num = str2double(last_training_path(strfind(last_training_path, 'ctl0')+4));
            end
        
        else 
            last_train_num = 0;
        end
        
        
        % identify probe type
        if contains(fpath, 'repeat')
            ltn = 7;
        elseif contains(fpath, 'four') | contains(fpath, 'two')
            ltn = 6;
        else
            error('check task type')
        end
        
        % find where to insert pre/post probe designator
        str_split = strfind(stage, 'probe')-1;
        
        
        if last_train_num == ltn %%%%% CHANGE
            stage = [stage(1:str_split) 'post' stage(str_split+1:end)];
        
            % number of days since last day of training
            numdays = datenum(current_train_date) - datenum(last_train_date);
            numdays_string = num2str(numdays);
            if length(numdays_string)<2
                numdays_string = ['0' numdays_string];
            elseif contains(numdays_string, '-')
                current_file = origin_paths_all{ifile}
                last_training_path
                error('negative numdays value')
            end
            
            % probe number
            switch 1
                case numdays >= 1 & numdays < 4
                    probe_num = '01';
                case numdays >= 4 & numdays < 10
                    probe_num = '02';
                case numdays >= 10 & numdays < 17
                    probe_num = '03';
                case numdays >= 17 & numdays < 24
                    probe_num = '04';
                case numdays >= 24 & numdays < 31
                    probe_num = '05';
                case numdays >= 31
                    probe_num = '06';
            end

        elseif last_train_num == 0
            stage = [stage(1:str_split) 'pre' stage(str_split+1:end)];
            probe_num = '01';
            numdays_string = ['01'];
        else
            
            % number of days since last day of training
            numdays = datenum(current_train_date) - datenum(last_train_date);
            numdays_string = num2str(numdays);
            if length(numdays_string)<2
                numdays_string = ['0' numdays_string];
            end
            
            stage = [stage(1:str_split) 'pre' stage(str_split+1:end)];
            probe_num = ['0' num2str(last_train_num+1)];
            
        end     

        save_file_name = [stage '_' probe_num '_' numdays_string 'd'];
        
        % skip if LED probe file and is already saved
        %if any(contains(get_file_paths_all(path_hold), save_file_name)) && contains(save_file_name, 'LED')
        %    continue
        %end
        if any(contains(get_file_paths_all(path_hold), save_file_name))
            continue
        end
        
    end
    
    %% reminder session
    if contains(save_file_name, 'rem')
        
        % info from current file
        [~, current_train_date, subject_id] = medass_mpc_id(origin_paths_all{ifile});
        
        % all origin paths containing subject id
        origin_paths_subj = origin_paths_all(contains(origin_paths_all, subject_id));
        
        % training days only       
        dest_paths_subj = get_file_paths_all(path_hold); % destination folder
        dest_paths_subj_train = dest_paths_subj(contains(dest_paths_subj, {'novar0', 'lovar0', 'mevar0', 'hivar0', 'exvar0', 'ctl0'}));
        origin_paths_subj_train = origin_paths_subj(contains(origin_paths_subj, {'novar0', 'lovar0', 'mevar0', 'hivar0', 'exvar0', 'ctl0'}));
        
        % iterate through destination ids finding corresponding origin
        % files
        comb_dest_origin_idx = false(size(origin_paths_subj_train,1),1);
        for idids = 1:size(dest_paths_subj_train,1)
            
            % training id
            dest_ids_1 = [dest_paths_subj_train{idids}(strfind(dest_paths_subj_train{idids},'novar0'):strfind(dest_paths_subj_train{idids},'novar0')+6) ...
                dest_paths_subj_train{idids}(strfind(dest_paths_subj_train{idids},'lovar0'):strfind(dest_paths_subj_train{idids},'lovar0')+6) ...
                dest_paths_subj_train{idids}(strfind(dest_paths_subj_train{idids},'mevar0'):strfind(dest_paths_subj_train{idids},'mevar0')+6) ...
                dest_paths_subj_train{idids}(strfind(dest_paths_subj_train{idids},'hivar0'):strfind(dest_paths_subj_train{idids},'hivar0')+6) ...
                dest_paths_subj_train{idids}(strfind(dest_paths_subj_train{idids},'exvar0'):strfind(dest_paths_subj_train{idids},'exvar0')+6) ...
                dest_paths_subj_train{idids}(strfind(dest_paths_subj_train{idids},'ctl0'):strfind(dest_paths_subj_train{idids},'ctl0')+4)];
            % session id
            dest_ids_2 = dest_paths_subj_train{idids}(strfind(dest_paths_subj_train{idids},'-')+1:strfind(dest_paths_subj_train{idids},'-')+2);
            
            % find corresponding origin files            
            comb_dest_origin_idx = comb_dest_origin_idx | contains(origin_paths_subj_train, [dest_ids_1 '\' dest_ids_2]) | contains(origin_paths_subj_train, [dest_ids_1 '_4t\' dest_ids_2])  | contains(origin_paths_subj_train, [dest_ids_1 '_2t\' dest_ids_2])...
                | contains(origin_paths_subj_train, [dest_ids_1 '_opto\' dest_ids_2]);
        end
        origin_paths_subj_train_indest = origin_paths_subj_train(comb_dest_origin_idx);

        % 
        if size(dest_paths_subj_train,1) ~= size(origin_paths_subj_train_indest,1)           
            origin_paths_subj_train_indest
            dest_paths_subj_train
            error('origin and destination file list mismatch')
        end
        
        % only include data files collected before probe
        path_dates = cell(size(origin_paths_subj_train_indest,1),1);
        for ipd = 1:size(origin_paths_subj_train_indest,1)
            [~,path_dates{ipd}] = medass_mpc_id(origin_paths_subj_train_indest{ipd});
            fclose('all');
        end        
        
        if ~isempty(path_dates) && ~isempty(dest_paths_subj_train(datenum(path_dates)<datenum(current_train_date)))
            
            % only dates before current session
            dest_paths_subj_train = dest_paths_subj_train(datenum(path_dates)<datenum(current_train_date));
            origin_paths_subj_train_indest = origin_paths_subj_train_indest(datenum(path_dates)<datenum(current_train_date));
            
            % search origin folder using info from destination folder
            last_training_path = origin_paths_subj_train_indest{end};
            [~, last_train_date] = medass_mpc_id(last_training_path);

            if contains(last_training_path, 'novar0')
                last_train_num = str2double(last_training_path(strfind(last_training_path, 'novar0')+6));
            elseif  contains(last_training_path, 'lovar0')
                last_train_num = str2double(last_training_path(strfind(last_training_path, 'lovar0')+6));
            elseif  contains(last_training_path, 'mevar0')
                last_train_num = str2double(last_training_path(strfind(last_training_path, 'mevar0')+6));
            elseif  contains(last_training_path, 'hivar0')
                last_train_num = str2double(last_training_path(strfind(last_training_path, 'hivar0')+6));
            elseif  contains(last_training_path, 'exvar0')
                last_train_num = str2double(last_training_path(strfind(last_training_path, 'exvar0')+6));
            elseif  contains(last_training_path, 'ctl0')
                last_train_num = str2double(last_training_path(strfind(last_training_path, 'ctl0')+4));
            end
        
        else 
            last_train_num = 0;
        end
        
        % number of days since last day of training
        numdays = datenum(current_train_date) - datenum(last_train_date);
        numdays_string = num2str(numdays);
        if length(numdays_string)<2
            numdays_string = ['0' numdays_string];
        elseif contains(numdays_string, '-')
            current_file = origin_paths_all{ifile}
            last_training_path
            error('negative numdays value')
        end

        % rem number
        rem_num = save_file_name(strfind(save_file_name,'-')+1:end);
        rem_num_str = num2str(rem_num);
        if length(rem_num_str)<2
            rem_num_str = ['0' rem_num_str];
        end

        save_file_name = [stage '_' rem_num_str '_' numdays_string 'd'];

             
        % skip if rem file is already saved
        if any(contains(get_file_paths_all(path_hold), save_file_name)) 
            continue
        end
        
    end
    
    
    
    

    %skip non-probe file if already saved
    %
    if exist([path_hold '\' save_file_name '.mat'] , 'file') && ~contains([path_hold '\' save_file_name '.mat'], 'probe')
        continue
    end
    
    
    % update name if needed
    %save_file_name = edit_paths(save_file_name);
    
    %load and save file
    load_medass_save(origin_paths_all{ifile}, path_hold, save_file_name);
    
    %close folders
    fclose('all');

end



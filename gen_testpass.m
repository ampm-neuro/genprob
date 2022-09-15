function [subject_ids, program_files, pass_fail, thirds_test_result] = gen_testpass(folder_path)
% loads med associates data files and determines whether the animal should
% be promoted to the next training stage
%
% input folder_path is the path of the folder containing the med associates
% data files to be tested. e.g., 'C:\Users\ampm\Desktop\ampm_testfolder' 
%
% gen_testpass outputs three variables. subject_ids is a cell array of 
% strings indicating the subject (mouse) id number of the subject that was
% run during that session. program_files is a cell array of strings
% indicating the med associates program file run during that session.
% pass_fail is a vector of 0's and 1's, whereby 0 indicates a fail and 1
% indicates a pass. 1's for promotion
%


%% pass criteria
min_np_rwds = 50;
min_train_trials = 42;
prop_rwd_gen_trials = [0.9 0.8 0.7 0.6 0.5 0.5];

%% get all files in folder path
file_paths = get_file_paths_all(folder_path);

% preallocate pass_fail
subject_ids = cell(size(file_paths,1),1);
program_files = cell(size(file_paths,1),1);
pass_fail = nan(size(file_paths,1),1);
thirds_test_result = nan(size(file_paths,1),1);

%% iterate through each file
for ipath = 1:size(file_paths,1)
    
   try % catch bad files

    % identify session type
        %file_paths{ipath}
        [program_file, ~, subject_id] = medass_mpc_id(file_paths{ipath});

        % load subject and program
        subject_ids{ipath} = subject_id;
        program_files{ipath} = program_file;

        % test each session type
        switch 1
        
            % if nose poke training
            case contains(program_file, 'NP')

                 % obtain session performance variables
                 % from med associates text file
                 fid = fopen(file_paths{ipath});
                 medass_output = textscan(fid, '%s%s%s%s%s%s');
                 for imav = 1:size(medass_output{1,1},1)

                     % check for counter variable information
                     if contains('A', medass_output{1,1}{imav,1}(1))...
                        && strcmp(medass_output{1,1}{imav,1}(2), ':')

                        %gather data related to that variable
                        perform_var_data = [];
                        for hvdi = 2:size(medass_output,2)                       
                            perform_var_data = [perform_var_data ...
                                medass_output{1,hvdi}(imav:imav+6)];
                        end
                        perform_var_data = perform_var_data';
                        perform_var_data = str2double(perform_var_data(:));
                        perform_var_data = perform_var_data(~isnan(perform_var_data));
                     end
                 end

                 % number of rewards
                 reward_ct = perform_var_data(12); 

                 % reward  minimum (defined above)
                 min_rwd = min_np_rwds;

                 % compare performance to required minimums
                 if reward_ct >= min_rwd
                     pass_fail(ipath) = 1;
                 elseif contains(program_file, 'NP01')
                     pass_fail(ipath) = 1;
                 else
                     pass_fail(ipath) = 0;
                 end
              
             %LED pretraining
             case contains(program_file, 'pretrain')

                % load file
                [medass_cell, trl_mtx] = load_medass(file_paths{ipath});
                
                % compute session performance measures
                trial_ct = size(trl_mtx,1);
                reward_ct = sum(~isnan(trl_mtx(:,11)));

                % performance minimums (partially defined above)
                min_trial = min_train_trials;
                min_rwd = trial_ct * 0.6 * prop_rwd_gen_trials(end);

                % compare performance to required minimums
                if trial_ct >= min_trial
                    pass_fail(ipath) = 1;
                else
                    pass_fail(ipath) = 0;
                end
                 

            % if gen training
            case contains(program_file, {'train', 'ctl', 'hivar', 'lovar', 'novar', 'mevar', 'exvar'})

                % load file
                %try
                    %file_paths{ipath}
                
                    [medass_cell, trl_mtx] = load_medass(file_paths{ipath});
                %catch
                    %disp(['error_file: ' file_paths{ipath}])
                    %error
                %end

                % test for minimum number of trials
                if size(trl_mtx,1) < min_train_trials
                    pass_fail(ipath) = 0;
                    continue % CHANGE TO PLOT PARTIAL OPTO SESSIONS
                end

                % identify the stage of gen training
                genstage_num = str2num(program_file(strfind(program_file, 'gen')+3 : strfind(program_file, 'gen')+4)); 

                % if pretraining
                if genstage_num < 7

                    % compute session performance measures
                    trial_ct = size(trl_mtx,1);
                    reward_ct = sum(~isnan(trl_mtx(:,11)));

                    % performance minimums (partially defined above)
                    min_trial = min_train_trials;
                    min_rwd = trial_ct * 0.6 * prop_rwd_gen_trials(genstage_num);

                    % compare performance to required minimums
                    if trial_ct >= min_trial && reward_ct >= min_rwd
                        pass_fail(ipath) = 1;
                    else
                        pass_fail(ipath) = 0;
                    end

                % if training
                elseif genstage_num >= 7

                    % compute wait durations
                    [wait_durations_all, wda_freq] = wait_times_prep(trl_mtx,2);
                    
                    % compute reward probabilities
                    [prob_dist, pd_freq] = rwd_prob_by_freq(medass_cell);
                    
                    rich_freq = pd_freq(prob_dist>0.5);
                    poor_freqs = pd_freq(prob_dist<0.5);

                    
                    
                    if isempty(rich_freq)
                        continue
                    end

                    % test differences between rich and poor tone wait times
                    %

                        [~, ~, anova_grp] = unique(wda_freq);
                        anova_p = anovan(wait_durations_all, anova_grp, 'display', 'off');
                                                

                        % posthoc ttests
                        tstat = [];
                        pvals = [];
                        for ipf = 1:length(poor_freqs)
                            [~, pval, ~, stats] = ttest2(wait_durations_all(wda_freq==poor_freqs(ipf)), wait_durations_all(wda_freq==rich_freq));
                            tstat = [tstat stats.tstat];
                            pvals = [pvals pval];
                        end
                        
                        % test
                        if all(tstat<0) && all(pvals<0.01) && anova_p<0.01 && length(wait_durations_all(wda_freq==rich_freq))>=2
                            pass_fail(ipath) = 1;
                        else
                            pass_fail(ipath) = 0;
                        end
                else
                    pass_fail(ipath) = nan;
                    warning('weird genstage_num')
                end

            %  probe session 
            case contains(program_file, 'probe')

                % load file
                [medass_cell, trl_mtx] = load_medass(file_paths{ipath});

                if size(trl_mtx,1) < min_train_trials
                    pass_fail(ipath) = 0;
                else
                    % check for nonflat preprobe
                    % (pvalue for nonflatness)
                    thirds_test_result(ipath) = thirds_test(trl_mtx);
                    %}
                    
                    pass_fail(ipath) = 1;
                end
                
                % if opto probe session
                if contains(program_file, 'opto')
                    if sum(trl_mtx,13) < 21
                        disp('insufficient opto trials')
                        pass_fail(ipath) = 0;
                    end
                end

           
            
            % if ambiguous data type   
            otherwise

                warning('No pass-fail test exists for this file type.')
                pass_fail(ipath) = 0;

        end
    
    catch
       file_paths{ipath}
    end
    
    % no figure for reminder sessions
    if contains(program_file, 'remind')
        continue
    end
        
    %figure
    %try

    figure;
        
            if contains(program_file, 'opto')
                                
                [~, trl_mtx] = load_medass(file_paths{ipath});
                
                wait_times_plot_opto(trl_mtx,3);
                
                % title 
                lastslash = strfind(file_paths{ipath}, 'Subject');
                pathstr = file_paths{ipath}(lastslash(end)+1+length('Subject'):end);
                pf_title = program_file;
                pf_title(strfind(pf_title,'_')) = [];
                pf_title(strfind(pf_title,'ampm'): strfind(pf_title,'ampm')+3) = [];
                title([pathstr ' ' pf_title ' ' num2str(pass_fail(ipath))] )
                ylim([0 60])
                
            elseif contains(program_file, 'probe')
                [~, trl_mtx] = load_medass(file_paths{ipath});
                wait_times_plot(trl_mtx,3);
                
                % title 
                lastslash = strfind(file_paths{ipath}, 'Subject');
                pathstr = file_paths{ipath}(lastslash(end)+1+length('Subject'):end);
                pf_title = program_file;
                pf_title(strfind(pf_title,'_')) = [];
                pf_title(strfind(pf_title,'ampm'): strfind(pf_title,'ampm')+3) = [];
                title([pathstr ' ' pf_title ' ' num2str(pass_fail(ipath))] )
                ylim([0 60])
            
            else
                
                try
                    
                [~, trl_mtx] = load_medass(file_paths{ipath});
                wait_times_plot(trl_mtx,3);
                
                % title 
                lastslash = strfind(file_paths{ipath}, 'Subject');
                pathstr = file_paths{ipath}(lastslash(end)+1+length('Subject'):end);
                pf_title = program_file;
                pf_title(strfind(pf_title,'_')) = [];
                pf_title(strfind(pf_title,'ampm'): strfind(pf_title,'ampm')+3) = [];
                title([pathstr ' ' pf_title ' ' num2str(pass_fail(ipath)) ' ' num2str([anova_p pvals])] )
                ylim([0 60])
                
                catch
                    close
                end
                
                
            end
            
            
    %catch
    %    close
    %end
    
    
end

%% Sort outputs

subject_ids

[subject_ids, sort_idx] = sort(subject_ids);
program_files = program_files(sort_idx);
pass_fail = pass_fail(sort_idx);
thirds_test_result = thirds_test_result(sort_idx,:);
pf = cell(size(pass_fail,1),1);
ttr = cell(size(thirds_test_result,1),1);
for ittr = 1:size(thirds_test_result,1)
    pf{ittr} = num2str(pass_fail(ittr));
    ttr{ittr} = num2str(thirds_test_result(ittr));
end
thirds_test_result = ttr;

[subject_ids program_files pf thirds_test_result]

drawnow
end


%% INTERNAL FUNCTIONS
%
function [mpc_id, start_date, subject_id] = medass_mpc_id(mpc_path)
%returns medass mpc file that was running when the medass data file at
%mpc_path was generated


%read output (table)
fid = fopen(mpc_path);
medass_output = textscan(fid, '%s%s%s%s%s%s');

%isolate session ID information
end_id = find(strcmp(medass_output{1,1}, 'MSN:'), 1, 'first');
mpc_id = [];
for mi = 1:size(medass_output,2)
    mpc_id = [mpc_id medass_output{1,mi}(end_id)];
end
mpc_id = mpc_id{find(contains(mpc_id, ':'))+1};


%isolate session start date information
end_id = find(strcmp(medass_output{1,1}, 'Start'), 1, 'first');
start_date = [];
for mi = 1:size(medass_output,2)
    start_date = [start_date medass_output{1,mi}(end_id)];
end
start_date = start_date{find(contains(start_date, ':'))+1};


%isolate subject id information
end_id = find(strcmp(medass_output{1,1}, 'Subject:'), 1, 'first');
subject_id = [];
for mi = 1:size(medass_output,2)
    subject_id = [subject_id medass_output{1,mi}(end_id)];
end
subject_id = subject_id{find(contains(subject_id, ':'))+1};

end

function [wait_times, frequencies] = wait_times_prep(trl_mtx, mean_or_all)
% outputs wait times and corresponding tone frequencies
% for means mean_or_all == 1, all == 2

probe_trials_idx = trl_mtx(:,3)==0;
all_wait_times = trl_mtx(probe_trials_idx,12)+2; %2s for fixed delay
all_tones = floor(trl_mtx(probe_trials_idx,2));

if mean_or_all == 1 % means only
    frequencies = unique(all_tones);
    wait_times = nan(size(frequencies));
    for itone = 1:length(frequencies)
        tone = frequencies(itone);
        wait_times(itone) = mean(all_wait_times(all_tones==tone));
    end
else % all waits
    wait_times = all_wait_times;
    frequencies = all_tones;
end
end

function item_paths = get_file_paths_all(folderpath)
%iterates through each folder in path and returns path incl folder

item_paths = cell(1);

%items in path
file_list = dir(folderpath);
for iitem = 1:length(file_list)
    current_sesh = file_list(iitem).name;
    
    %omited folder and file names (folders named 'old' are invisible)
    if strcmp(current_sesh, '.') || strcmp(current_sesh, '..') || strcmp(current_sesh, 'old')
        continue
    end

    if isfolder([folderpath '\' file_list(iitem).name])
        item_paths_hold = get_file_paths_all([folderpath '\' file_list(iitem).name]);
        item_paths = [item_paths; item_paths_hold];
    elseif isfile([folderpath '\' file_list(iitem).name])
        item_paths = [item_paths; {[folderpath '\' file_list(iitem).name]}];
    end

end

%remove empty cells
item_paths = item_paths(find(~cellfun(@isempty, item_paths)));
end

function [rwd_prob, unq_frq] = rwd_prob_by_freq(medass_cell)
% outputs the probability of reward for each frequency

% raw frequencies
freq_raw = floor(medass_cell{21});

% raw reward probablities
rwd_probs_raw = medass_cell{22};

% unique freq
unq_frq = unique(freq_raw);

% compute reward probabilities
rwd_prob = nan(size(unq_frq));
for ifrq = 1:length(unq_frq)
    rwd_prob(ifrq) = mean(rwd_probs_raw(freq_raw==unq_frq(ifrq)));
end

% ensure proportion (not probability)
if rwd_prob(1)>1
    rwd_prob = rwd_prob./100;
end
end



function trl_time_out = trl_time(cell_of_interest, trl_start, trl_end)
    trl_time_out = cell_of_interest>=trl_start & cell_of_interest<=trl_end;
end





function probe_problem_interact_plots(subject_id, training_condition)
% a series of plot showing how problem learning is affected by and affects
% probe responses

% training_condition = 'train_mevar';

% subject folder paths
overall_folder_name = 'C:\Users\ampm1\Documents\MATLAB\generalization\behavioral_data\two_tone\';
folder_name = [overall_folder_name training_condition];
subject_paths = get_file_paths_targeted(folder_name, subject_id);

% plots
%
    ProbeProbSesh_comp = [1 1 1; 2 1 2; 2 2 1; 3 2 2 ; 3 3 1; 4 3 2; 4 4 1 ; 5 4 2; 5 5 1; 6 5 2; 6 6 1; 7 6 2]; 

for iprobe = 7:-1:1
    for iproblem = 6:-1:1
        for isesh = 2:-1:1
           
            if ~any(ismember(ProbeProbSesh_comp, [iprobe iproblem isesh], 'rows'))
                continue
            end

            
            % probe
            if iprobe < 7
                probe_path = subject_paths(contains(subject_paths, ['preprobe_0' num2str(iprobe)]));
            elseif iprobe >= 7
                probe_path = subject_paths(contains(subject_paths, ['postprobe_0' num2str(iprobe-6)]));
            end
            
            
            % check if probe completed
            if isempty(probe_path)
                continue
            else
                figure; hold on
            end
            
            probe_path = probe_path{1};
            load(probe_path, 'trl_mtx')
            wait_times_plot(trl_mtx,1);
            
            
            % problem
            problem_path = subject_paths(contains(subject_paths, ['var0' num2str(iproblem)]) | contains(subject_paths, ['ctl0' num2str(iproblem)]));
            if isesh==1
                problem_path = problem_path{1};
            elseif isesh==2
                problem_path = problem_path{end};
            end
            load(problem_path, 'trl_mtx')
            wait_times_plot(trl_mtx,3, [0 0 0]);
            
            ylim([0 40])
            if isesh==1
                title(['probe ' num2str(iprobe) ', first session of problem ' num2str(iproblem)])
            else
                title(['probe ' num2str(iprobe) ', last session of problem ' num2str(iproblem)])
            end
        end
    end
end



function [all_fixed_effects, all_beta_pvals, all_wait_times] = noplot_allprobes(file_paths, varargin)
% plot all preprobes

%input
if nargin==2
    interp_option = varargin{1};
else
    interp_option = 0;
end

if isempty(file_paths)
    all_fixed_effects = [];
    all_beta_pvals = [];
    all_wait_times = [];
    return
end

% identify probe numbers
pp_nums = nan(size(file_paths,1),1);
for ipath = 1:size(file_paths,1)
   
    u_idx = strfind(file_paths{ipath}, '_');
    probe_num_string = file_paths{ipath}(u_idx(end-1)+1 : u_idx(end)-1);
    
    pp_nums(ipath) = str2double(probe_num_string);
    
end

% preallocate
all_fixed_effects = [];
all_beta_pvals = [];
all_wait_times = [];

% iterate through unique probe numbers
legend_strings = cell(1, max(pp_nums));
for inum = 1:max(pp_nums)
    
    ppnum_fps = file_paths(pp_nums==inum); 
    
    %plot curves
    %{
    try
        if isempty(ppnum_fps)
            plot([1 1+realmin],[1 1], 'color', colors(inum,:))
        elseif size(ppnum_fps,1)==1
            load(ppnum_fps{1}, 'trl_mtx')
            plot_normal_fit_subj(trl_mtx, 2, colors(inum,:));
            clear trl_mtx
        else
            [fixed_effects, random_effects, beta_pvals, coef_names, modelFun, unqfrq] = ...
                mixed_model_fit(ppnum_fps, 1, colors(inum,:));
        end
        
        all_fixed_effects = [all_fixed_effects; fixed_effects'];
        all_beta_pvals = [all_beta_pvals; beta_pvals];
        
        
        catch
    end
    
    %}
    
    
    %plot means
    %
    if ~isempty(ppnum_fps)
        
        %get wait times from every session
        wait_times_all = [];
        frequencies_all = [];
        for ippnum = 1:size(ppnum_fps,1)
            load(ppnum_fps{ippnum}, 'trl_mtx')
            
            try
                [wait_times_local, frequencies_local] = wait_times_prep(trl_mtx, 1, interp_option);
            catch
                disp('used old')
                [wait_times_local, frequencies_local] = wait_times_prep(trl_mtx, 1, 1);
            end
            
            
            if size(wait_times_local,2)>size(wait_times_local,1)
                wait_times_local = wait_times_local';
            end
            
            %old = [size(wait_times_local); size(frequencies_local)]
            
            
            %new = [size(wait_times_local); size(frequencies_local)]
            
            
            %wait_times_local = zscore_mtx(wait_times_local);
            wait_times_all = [wait_times_all; wait_times_local];
            frequencies_all = [frequencies_all; frequencies_local];
            all_wait_times = [all_wait_times; wait_times_local'];
        end
        
        % means and ses
        unq_frq = unique(frequencies_all);
        mean_wait_times = nan(length(unq_frq),1);
        se_wait_times = nan(size(mean_wait_times));
        for ifrq = 1:length(unq_frq)
            mean_wait_times(ifrq) = nanmean(wait_times_all(frequencies_all==unq_frq(ifrq)));
            se_wait_times(ifrq) = nanstd(wait_times_all(frequencies_all==unq_frq(ifrq)))/sqrt(sum(frequencies_all==unq_frq(ifrq)));
        end

    else
        
    end
    %}
    
    
end


function [wait_times, frequencies, freq_numbers, modelFun, coefEsts] = wait_times_prep(trl_mtx, mean_or_all, varargin)
% [wait_times, frequencies] = wait_times_prep(trl_mtx, mean_or_all, interp_input)
% outputs wait times and corresponding tone frequencies
% for means mean_or_all == 1, all == 2
% interp_input defaults to off, input 1 to guess missing values



load('unqfrq41.mat', 'unqfrq41')

if nargin == 3
    interp_input = varargin{1};
else
    interp_input = -1;
end

probe_trials_idx = trl_mtx(:,3)==0;
all_wait_times = trl_mtx(probe_trials_idx,12)+2; %2s for fixed delay
all_tones = trl_mtx(probe_trials_idx,2);
all_tones = floor(all_tones);
modelFun = [];
coefEsts = [];

% remove wait times greater than 60s, since the tone stopped on those
% trials
max_wait_time = 60; %s
max_wait_elim = all_wait_times>max_wait_time;
all_wait_times(max_wait_elim) = [];
all_tones(max_wait_elim) = [];


% zscore waits
%all_wait_times = zscore_mtx(all_wait_times);

if mean_or_all == 1 % means onlynbnb
    frequencies = unique(all_tones);
    wait_times = nan(size(frequencies));
    for itone = 1:length(frequencies)
        tone = frequencies(itone);
        wait_times(itone) = mean(all_wait_times(all_tones==tone));
    end
else % all waits
    wait_times = all_wait_times;
    frequencies = all_tones;
    
    % sort by frequency
    [frequencies, sort_idx] = sort(frequencies);
    wait_times = wait_times(sort_idx);
end


% model-fit interpolation and extrapolation
if interp_input == 0
    all_freq_nums = 1:length(unqfrq41);
    obs_freq_nums = all_freq_nums(ismember(unqfrq41, unique(all_tones))); 
    try
    [~, coefEsts, modelFun] = ampm_normal_logistic_fit(obs_freq_nums, wait_times, [nanmean(wait_times) 21 0 1]);
    catch
        try
        [~, coefEsts, modelFun] = ampm_normal_logistic_fit(obs_freq_nums, wait_times, [nanmean(wait_times) 20.5 0 1]);
        catch
            try
            [~, coefEsts, modelFun] = ampm_normal_logistic_fit(obs_freq_nums, wait_times, [nanmean(wait_times) 21.5 0 1]);
            catch
                try
                [~, coefEsts, modelFun] = ampm_normal_logistic_fit(obs_freq_nums, wait_times, [nanmean(wait_times) 20.0 0 1]);
                catch
                    try
                    [~, coefEsts, modelFun] = ampm_normal_logistic_fit(obs_freq_nums, wait_times, [nanmean(wait_times) 22.0 0 1]);
                    catch
                    end
                end
            end
        end
    end
    wait_times_hold = nan(1, length(all_freq_nums));
    wait_times_hold(ismember(all_freq_nums, obs_freq_nums)) = wait_times;
    wait_times_hold(~ismember(all_freq_nums, obs_freq_nums)) = modelFun(coefEsts, setdiff(all_freq_nums, obs_freq_nums));
    wait_times = wait_times_hold;
    frequencies = unqfrq41';

% interpolation
elseif interp_input==1
    all_freq_nums = 1:length(unqfrq41);
    obs_freq_nums = all_freq_nums(ismember(unqfrq41, unique(all_tones)));
    
    waits_interp = interp1(obs_freq_nums, wait_times, all_freq_nums, 'linear');
    
    wait_times = waits_interp';
    wait_times(wait_times<2.01) = 2.01;
    frequencies = unqfrq41';


% interpolation & extrapolation
elseif interp_input==2
    all_freq_nums = 1:length(unqfrq41);
    obs_freq_nums = all_freq_nums(ismember(unqfrq41, unique(all_tones)));

    wait_times_input = nan(size(all_freq_nums));
    wait_times_input(obs_freq_nums) = wait_times;
    
    [wait_times_smooth] = nansmooth_adaptive_ampm(wait_times_input, 7, 5);
    wait_times_smooth(obs_freq_nums) = wait_times;
    wait_times = wait_times_smooth';

    wait_times(wait_times<2.01) = 2.01;
    frequencies = unqfrq41';
    

% interpolation & extrapolation & smooth
elseif interp_input==3
    all_freq_nums = 1:length(unqfrq41);
    obs_freq_nums = all_freq_nums(ismember(unqfrq41, unique(all_tones)));

    wait_times_input = nan(size(all_freq_nums));
    wait_times_input(obs_freq_nums) = wait_times;
    
    [wait_times_smooth] = nansmooth_adaptive_ampm(wait_times_input, 3, 3);
    wait_times = wait_times_smooth';

    wait_times(wait_times<2.01) = 2.01;
    frequencies = unqfrq41';
    
    
% raw with nans
elseif interp_input==4
    
    all_freq_nums = 1:length(unqfrq41);
    obs_freq_nums = all_freq_nums(ismember(unqfrq41, unique(all_tones)));
    wt_hold = wait_times;
    wait_times = nan(size(unqfrq41))';
    wait_times(obs_freq_nums) = wt_hold;
    frequencies = unqfrq41';
    
end

% frequency numbers
freq_numbers = nan(length(frequencies),1); 
for ifrq = 1:length(freq_numbers)
    freq_numbers(ifrq) = find(ismember(unqfrq41,frequencies(ifrq)));
end




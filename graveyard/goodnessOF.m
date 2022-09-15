function gofz = goodnessOF(obs_data, fit_data)
% Goodness of fit between test and reference data. 
% gof = goodnessof(obs_data, fit_data, cost_func) returns the goodness of 
% fit between the data, obs_data, and the reference, fit_data.
%
% Normalized root-mean-square deviation
if any(size(obs_data) ~= size(fit_data))
    error('obs_data and fit_data must be same size')
end

%remove nans
nnan_idx = ~isnan(obs_data) & ~isnan(fit_data);
obs_data = obs_data(nnan_idx);
fit_data = fit_data(nnan_idx);

%gof
%gofz = immse(obs_data, fit_data) / immse(obs_data, repmat(mean(obs_data), size(obs_data)))
%gofz = immse(obs_data, fit_data) / variance(obs_data)
%gofz = stand_dev(obs_data - fit_data) / mean(obs_data)

% std relative to model divided by std relative to mean
% Normalized root-mean-square deviation
% no standard way to normalize... many use range or mean
gofz = stand_dev(obs_data - fit_data) ./ stand_dev(obs_data - repmat(mean(obs_data), size(obs_data)));


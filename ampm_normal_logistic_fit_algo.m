function [fit_yvals, coeffs_out, modelFun_full] = ampm_normal_logistic_fit_algo(xpos, ypos, varargin)
% computes the mean, std, multiplier, and intercept of the best fit normal 
% distribution

% model
% p(1) = intercept
% p(2) = inflection tonal position (norm and logistic)
% p(3) = logistic multiplier (STM strength)
% p(4) = normal multiplier (LTM strength)

% acceptable bounds
%acceptable_coef_bounds = [-5 45; -5 41; -20 20; -20 20];
%acceptable_coef_bounds = [-30 60; -30 41; -30 30; -30 30];
acceptable_coef_bounds = [-30 60; -30 71; -30 30; -30 30];

%% common model details 

% starting vals

    if nargin == 3
        startingVals = varargin{1};
    else
        startingVals = [nanmean(ypos) 20 0 1];
    end
    
% alpha level for determining if one model is better than another (try
% 0.50)
%alpha_level = 0.05;
alpha_level = 0.999;


%% model specific details

% intercept only
    
    % p(1) = intercept
    startingVals_int = startingVals(1);
    modelFun_int = @(p,x) p(1);

% intercept and logistic

    % model
    % p(1) = intercept
    % p(2) = inflection tonal position (norm and logistic)
    % p(3) = logistic multiplier (STM strength)
    startingVals_log = startingVals([1 2 3]);
    modelFun_log = @(p,x) (p(3).*(exp(x-p(2))./(exp(x-p(2))+1))) + p(1);


% intercept and normal

    % model
    % p(1) = intercept
    % p(2) = inflection tonal position (norm and logistic)
    % p(3) = normal multiplier (LTM strength)
    startingVals_norm = startingVals([1 2 4]);
    modelFun_norm = @(p,x) p(1) + p(3).*exp(-((x-p(2))./abs(p(3))^(3/4)).^2);


% intercept, logistic, and normal

    % model
    % p(1) = intercept
    % p(2) = inflection tonal position (norm and logistic)
    % p(3) = logistic multiplier (STM strength)
    % p(4) = normal multiplier (LTM strength)
    startingVals_full = startingVals;
    modelFun_full = @(p,x) (p(4).*exp(-((x-p(2))./(abs(p(4)))^(3/4)).^2))...
        + (p(3).*(exp(x-p(2))./(exp(x-p(2))+1))) + p(1);
    
    % modelFun_full = @(p,x) (p(4).*exp(-((x-p(2))./(abs(p(4)))^(1/2)).^2))...
    %    + (p(3).*(exp(x-p(2))./(exp(x-p(2))+1))) + p(1);


    

%% fit each possible model; compute MSE vectors

% make uniform vectors
nan_hold = nan(1,41); 
nan_hold(xpos) = ypos; 
ypos = nan_hold;
xpos = 1:41;


% intercept only
MSE_int = nan(size(ypos));
for idp = 1:length(ypos)
    MSE_int(idp) = (ypos(idp)-nanmean(ypos(setdiff(1:length(ypos),idp)))).^2;
end

% intercept and logistic
options_log = statset('MaxIter', 1000000,'FunValCheck','on');
MSE_log = whaf(xpos, ypos, modelFun_log, startingVals_log, startingVals, acceptable_coef_bounds, options_log);

% intercept and normal
options_norm = statset('MaxIter', 1000000,'FunValCheck','on');
MSE_norm = whaf(xpos, ypos, modelFun_norm, startingVals_norm, startingVals, acceptable_coef_bounds, options_norm);

% intercept, logistic, and normal
options_full = statset('MaxIter', 1000000,'FunValCheck','on');
MSE_full = whaf(xpos, ypos, modelFun_full, startingVals_full, startingVals, acceptable_coef_bounds, options_full);


%[nanmean(MSE_int) nanmean(MSE_log) nanmean(MSE_norm) nanmean(MSE_full)]


%% compare models (algorithmically)
%figure; errorbar_barplot([{(MSE_int)} {(MSE_log)} {(MSE_norm)} {(MSE_full)}])

% keep one of normal or logistic model with smaller mean MSE
    if nanmean(MSE_log) <= nanmean(MSE_norm)
        %disp('log better than norm')

        % if the logistic model is better fit than the intercept model
        if ~isnan(ttest(MSE_int, MSE_log)) && ttest(MSE_int, MSE_log, 'alpha', alpha_level) && nanmean(MSE_log) <= nanmean(MSE_int)
            %disp('log better than int')

            % if full model is better fit than the logistic model
            if ~isnan(ttest(MSE_log, MSE_full)) && ttest(MSE_log, MSE_full, 'alpha', alpha_level) && nanmean(MSE_full) <= nanmean(MSE_log)
                %disp('full better than log')

                % full model
                coeffs_out = nlinfit(xpos, ypos, modelFun_full, startingVals_full, options_full);

            else
                %disp('log better than full')
                
                % logistic only
                coeffs_out = nlinfit(xpos, ypos, modelFun_log, startingVals_log, options_log);
                coeffs_out = [coeffs_out(1:3) 0];
            end

        else
            %disp('int better than log')
            % intercept only
            coeffs_out = [startingVals(1) 0 0 0];
        end

    else

        % if the normal model is better fit than the intercept model
        if ~isnan(ttest(MSE_int, MSE_norm)) && ttest(MSE_int, MSE_norm, 'alpha', alpha_level) && nanmean(MSE_norm) <= nanmean(MSE_int)

            % if full model is better fit than the logistic model
            if ~isnan(ttest(MSE_norm, MSE_full)) && ttest(MSE_norm, MSE_full, 'alpha', alpha_level) && nanmean(MSE_full) <= nanmean(MSE_norm)

                % full model
                coeffs_out = nlinfit(xpos, ypos, modelFun_full, startingVals_full, options_full);

            else
                % normal only
                coeffs_out = nlinfit(xpos, ypos, modelFun_norm, startingVals_norm, options_norm);
                coeffs_out = [coeffs_out(1:2) 0 coeffs_out(3)];
            end

        else
            % intercept only
            coeffs_out = [startingVals(1) 0 0 0];
        end

    end




%% compute modeled yvalues
fit_yvals = modelFun_full(coeffs_out, unique(xpos));
  


end

%% withold and fit function
function MSE = whaf(xpos, ypos, modelFun, startingVals_model, startingVals_all, acceptable_coef_bounds, options)
% iteratively remove 1 data point, fit, compute MSE of witheld point from
% the fit model. Output vector of MSEs.

    % preallocate modeled values of witheld data
    wp_val = nan(size(ypos));
    
    % iterate through data points
    for idp = xpos

        if isnan(xpos)
            continue
        end
        
        % remaining data points
        rp = setdiff(1:length(xpos),idp);

        %try
            % fit model of remaining data points
            fm_coefs = nlinfit(xpos(rp), ypos(rp), modelFun, startingVals_model, options);
            
            test_coefs_in = zeros(1,4);
            test_coefs_in(ismember(startingVals_all, startingVals_model)) = fm_coefs;
            if test_coefs(test_coefs_in, acceptable_coef_bounds)
            
                % out of range coefs catch
                wp_val(idp) = nan;
            
            else
            
                % compute value of model at witheld point 
                wp_val(idp) = modelFun(fm_coefs, idp);
                
            end

        %catch
            % error catch
        %    display('ugh')
        %    wp_val(idp) = nan;
        %end

    end

    % compute MSE
    MSE = (ypos - wp_val).^2;
  
end


%% test coef bounds
function testpass = test_coefs(fm_coefs, coef_bounds)
    
    testpass = any(...
                    fm_coefs(1) < coef_bounds(1,1) | fm_coefs(1) > coef_bounds(1,2) | isnan(fm_coefs(1)) | ...
                    fm_coefs(2) < coef_bounds(2,1) | fm_coefs(2) > coef_bounds(2,2) | isnan(fm_coefs(2)) |...
                    fm_coefs(3) < coef_bounds(3,1) | fm_coefs(3) > coef_bounds(3,2) | isnan(fm_coefs(3)) |...
                    fm_coefs(4) < coef_bounds(4,1) | fm_coefs(4) > coef_bounds(4,2) | isnan(fm_coefs(4)) ...
                    );
end


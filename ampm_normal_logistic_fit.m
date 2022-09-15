function [fit_yvals, coefEsts, modelFun, residuals] = ampm_normal_logistic_fit(xpos, ypos, varargin)
% computes the mean, std, multiplier, and intercept of the best fit normal 
% distribution



% starting vals
if nargin == 3
    startingVals = varargin{1};
else
    startingVals = [nanmean(ypos) 21 0 1];
end

% model
% p(1) = intercept
% p(2) = logistic multiplier (short-term mem strength)
% p(3) = logistic x postion
% p(4) = normal multiplier (long-term mem strength)
% p(5) = mean (long-term mem accuracy)
% p(6) = std (long-term mem specificity)
%{
modelFun = @(p,x) (p(4).*exp(-((x-p(5))./p(6)).^2))...
    + (p(2).*(exp(x-p(3))./(exp(x-p(3))+1)))...
    + p(1);
%}
%
%
% model
% p(1) = intercept
% p(2) = inflection tonal position (norm and logistic)
% p(3) = logistic multiplier (STM strength)
% p(4) = normal multiplier (LTM strength)

%{
modelFun = @(p,x) (p(4).*exp(-((x-p(2))./(abs(p(4)))^(3/4)).^2))...
    + (p(3).*(exp(x-p(2))./(exp(x-p(2))+1)))...
    + p(1);

%}
modelFun = @(p,x) (p(4).*exp(-((x-p(2))./abs(p(4))^(3/4)).^2))...
    + (p(3).*(exp(x-p(2))./(exp(x-p(2))+1)))...
    + p(1);

% orientation
if size(xpos,2)>size(xpos,1)
    if size(ypos,2)<=size(ypos,1)
        ypos = ypos';
    end
elseif size(xpos,2)<size(xpos,1)
    if size(ypos,2)>=size(ypos,1)
        ypos = ypos';
    end
end


% fit
options = statset('MaxIter',10000000000,'FunValCheck','on');
[coefEsts, residuals] = nlinfit(xpos, ypos, modelFun, startingVals, options);
coefEsts = coefEsts;

% output
fit_yvals = modelFun(coefEsts, unique(xpos));


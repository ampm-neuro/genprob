function [fit_yvals, coefEsts, modelFun, GoF] = ampm_logistic_fit(xpos, ypos, varargin)
% computes the mean, std, multiplier, and intercept of the best fit normal 
% distribution

% model
% p(1) = intercept
% p(2) = inflection tonal position (norm and logistic)
% p(3) = logistic multiplier (STM strength)
%modelFun = @(p,x) (p(3).*(exp(x-p(2))./(exp(x-p(2))+1))) + p(1);
modelFun = @(p,x) p(3)./(1+exp(p(2)-x)) + p(1);

% fit
if nargin==3
    startingVals = varargin{1};
else
    startingVals = [nanmean(ypos) 20 0];
end
options = statset('MaxIter',1000000,'FunValCheck','on');
coefEsts = nlinfit(xpos, ypos, modelFun, startingVals, options);

% output
fit_yvals = modelFun(coefEsts, unique(xpos));
GoF = [];%goodnessOF(ypos, fit_yvals);

function [fit_yvals, coefEsts, modelFun, GoF] = ampm_normal_fit(xpos, ypos, startingVals)
% computes the mean, std, multiplier, and intercept of the best fit normal 
% distribution

% model
% p(1) = intercept
% p(2) = multi (mem strength)
% p(3) = mean (accuracy)
% p(4) = std (specificity)
modelFun = @(p,x) p(2).*exp(-((x-p(3))./p(4)).^2) + p(1);

% fit
%startingVals = [12.6 10 21 9];
%startingVals = [1 1 1 1].*20;
options = statset('MaxIter',1000000,'FunValCheck','on');
coefEsts = nlinfit(xpos, ypos, modelFun, startingVals, options);

% output
fit_yvals = modelFun(coefEsts, unique(xpos));
GoF = [];%goodnessOF(ypos, fit_yvals);

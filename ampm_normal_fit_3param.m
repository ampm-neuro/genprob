function [fit_yvals, coefEsts, modelFun, GoF] = ampm_normal_fit_3param(xpos, ypos, startingVals)
% computes the mean, std, multiplier, and intercept of the best fit normal 
% distribution

% model
% p(1) = intercept
% p(2) = height and std
% p(3) = mean (accuracy)
modelFun = @(p,x) p(1) + p(2).*exp(((x-p(3)).*(p(3)-x))./abs(p(2)).^(3/2));


% fit
startingVals = [12.6 10 21 9];
%startingVals = [1 1 1 1].*20;
options = statset('MaxIter',1000000,'FunValCheck','on');
coefEsts = nlinfit(xpos, ypos, modelFun, startingVals, options);

% output
fit_yvals = modelFun(coefEsts, unique(xpos));
GoF = [];%goodnessOF(ypos, fit_yvals);

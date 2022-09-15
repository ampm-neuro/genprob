function [fit_yvals, coefEsts, modelFun] = ampm_timecurve_fit(xpos, ypos, startingVals)
% computes the mean, std, multiplier, and intercept of the best fit normal 
% distribution

% normal model
% p(1) = intercept
% p(2) = multi
% p(3) = mean (accuracy)
% p(4) = std (specificity)
 modelFun = @(p,x) p(1) + p(2).*exp(-((x-p(3))./p(4)).^2);

% gamma model
% p(1) = intercept
% p(2) = multi
% p(3) = mean (accuracy)
% p(4) = std (specificity)
%modelFun = @(p,x) p(1)+ (p(2).*gampdf(x, p(3), p(4)));

% linear fit
%p(1) = intercept
%p(2) = slope
%modelFun = @(p,x) p(1)+ (p(2).*x) + x.^p(3)


% fit
%options = statset('MaxIter',1000000,'FunValCheck','on');
%coefEsts = nlinfit(xpos, ypos, modelFun, startingVals, options);
coefEsts = lsqcurvefit(modelFun,startingVals,xpos,ypos, 0);

% output
fit_yvals = modelFun(coefEsts, unique(xpos));

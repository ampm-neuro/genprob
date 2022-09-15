function [fit_yvals, coefEsts, GoF] = ampm_comb_fit(xpos, ypos)
% computes the mean, std, multiplier, and intercept of the best fit normal 
% distribution
xpos = xpos(:);
ypos = ypos(:);

% fit requires no nans
nnan_idx = ~isnan(xpos) & ~isnan(ypos);
xpos_nnan = xpos(nnan_idx);
ypos_nnan = ypos(nnan_idx);

% model options
%{
% old options
opts = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[1 0.01 0 -1 0],...
               'Upper',[41 41 100 1 100],...
               'Startpoint',[15 2 14 0 11],...
               'MaxIter', 100,...
               'TolFun', 0.00000001,...
               'TolX', 0.00000001);
%}
           
%{           
opts = fitoptions('Method','NonlinearLeastSquares',...
   'Lower',[1 0.05 1 -1 0],...
   'Upper',[41 41 100 1 100],...
   'Startpoint',[15 2 20 0 11],...
   'MaxIter', 100,...
   'TolFun', 0.000000001,...
   'TolX', 0.000000001);
           
           
           
ft = fittype('((exp(-0.5 * ((x - a)./b).^2) ./ (sqrt(2*pi) .* b) - min((exp(-0.5 * ((x - a)./b).^2) ./ (sqrt(2*pi) .* b))))./(max(exp(-0.5 * ((x - a)./b).^2) ./ (sqrt(2*pi) .* b)) - min(exp(-0.5 * ((x - a)./b).^2) ./ (sqrt(2*pi) .* b)))).*c + d.*x + e','options',opts);

% fit coefficients
FO = fit(xpos_nnan, ypos_nnan, ft);
coefEsts(1) = FO.a;
coefEsts(2) = FO.b;
coefEsts(3) = FO.c;
coefEsts(4) = FO.d;

% evaluate at original x positions
fit_yvals = feval(FO, xpos);

% goodness of fit
GoF = goodnessOF(ypos, fit_yvals);

coefEsts
%}
% stats
%modelFun = @(coefEsts,x) ((exp(-0.5 * ((x - coefEsts(1))./coefEsts(2)).^2) ./ (sqrt(2*pi) .* coefEsts(2)) - min((exp(-0.5 * ((x - coefEsts(1))./coefEsts(2)).^2) ./ (sqrt(2*pi) .* coefEsts(2)))))./(max(exp(-0.5 * ((x - coefEsts(1))./coefEsts(2)).^2) ./ (sqrt(2*pi) .* coefEsts(2))) - min(exp(-0.5 * ((x - coefEsts(1))./coefEsts(2)).^2) ./ (sqrt(2*pi) .* coefEsts(2))))).*coefEsts(3) + coefEsts(4).*x + coefEsts(5);
%fitnlm(xpos_nnan,ypos_nnan,modelFun,[15 2 20 0 11], 'CoefficientNames', {'mu', 'sigma', 'mult', 'slope', 'intercept'})

modelFun = @(coefEsts,x)  coefEsts(1) + coefEsts(2).*x +coefEsts(3).*((exp(-0.5 * ((x - coefEsts(4))./coefEsts(5)).^2) ./ (sqrt(2*pi) .* coefEsts(5)) - min((exp(-0.5 * ((x - coefEsts(4))./coefEsts(5)).^2) ./ (sqrt(2*pi) .* coefEsts(5)))))./(max(exp(-0.5 * ((x - coefEsts(4))./coefEsts(5)).^2) ./ (sqrt(2*pi) .* coefEsts(5))) - min(exp(-0.5 * ((x - coefEsts(4))./coefEsts(5)).^2) ./ (sqrt(2*pi) .* coefEsts(5)))));
FO = fitnlm(xpos_nnan,ypos_nnan,modelFun,[11 0 20 15 2], 'CoefficientNames', {'intercept', 'slope', 'mult', 'mu', 'sigma'});
coefEsts = FO.Coefficients.Estimate';

GoF = goodnessOF(ypos, feval(FO, xpos));

fit_yvals = feval(FO, xpos);


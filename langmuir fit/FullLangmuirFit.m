function [fitresult, gof] = FullLangmuirFit(t, filtered, concentration, scaleMax, stopLower, stopUpper)
%CREATEFIT(T,FILTERED)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : t
%      Y Output: filtered
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 21-Jul-2020 16:15:22


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( t, filtered );

% Set up fittype and options.
ft = fittype( 'a*(kon*c/(kon*c+koff)).*(1-exp(-(kon*c+koff).*(x))).*(x <= stop)+((kon*c/(kon*c+koff)) * (1 - exp(-(kon*c+koff) * stop)) * exp(koff * stop))*a*exp(-koff * (x)).*(x > stop);', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0.2 concentration 1e-6 1000 stopLower];
opts.StartPoint = [0.2 concentration 1e-5 1e5 stopLower];
opts.Upper = [scaleMax concentration 1e-4 1000000 stopUpper];

%opts.BOUND = [scale concentration koff kon stoptimeAssociation]

%lsqcurvefit

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
% 
% %Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'filtered vs. t', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel t
% ylabel filtered
% grid on



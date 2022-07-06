function [fitresult, gof] = createFit3(Xaxes, SNRmedian)
%CREATEFIT(XAXES,SNRMEDIAN)
%  Create a fit.
%
%  Data for 'powerFit' fit:
%      X Input : Xaxes
%      Y Output: SNRmedian
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 20-Jul-2020 12:08:51


%% Fit: 'powerFit'.
[xData, yData] = prepareCurveData( Xaxes, SNRmedian );

% Set up fittype and options.
ft = fittype( 'a*x^b', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0 0.210609553933221];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'powerFit' );
h = plot( fitresult, xData, yData );
legend( h, 'SNRmedian vs. Xaxes', 'powerFit', 'Location', 'NorthEast' );
% Label axes
xlabel Xaxes
ylabel SNRmedian
grid on



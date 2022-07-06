function y_est = langmuirModel(kon, koff, c, a, x, stop)


    y_est = a*(kon*c/(kon*c+koff)).*(1-exp(-(kon*c+koff).*(x))).*(x <= stop)+((kon*c/(kon*c+koff)) * (1 - exp(-(kon*c+koff) * stop)) * exp(koff * stop))*a*exp(-koff * (x)).*(x > stop);

end
function [curve,y_assoc,y_dissoc] = langmuirModel(k_on, k_off, smax, scale, timePoints, stopTime)
%     k_on = 2.28e4; 
%     k_off = 2.36e-05;
    concentration = 70e-9;
%     smax =  1;
%     scale =  100;

    t = timePoints;
    stop_time = stopTime;
    stop_time_ind = find(t == stop_time);
    %%

    beta = k_on * concentration * smax;
    alpha = k_on * concentration + k_off;

    y_0 = beta / alpha * (1 - exp(-alpha * stop_time)) * exp(k_off * stop_time);
    y_assoc = (scale*(beta / alpha).*(1 - exp(-alpha * t(1:stop_time_ind))));
    y_dissoc = scale*y_0.*exp(-k_off * t(stop_time_ind+1:end));

    curve = [y_assoc y_dissoc];
    %plot(t, curve)






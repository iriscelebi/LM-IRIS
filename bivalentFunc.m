

%def bivalent(self, timePoints, p0,p1,p2,p3,p4, concentration, smax, t_stop)
function newResult = bivalentFunc(k_on1,k_off1,k_on2,k_off2,B0, A0, scale, t_off,time)

%     timePoints = jamesData(:,1);


%     k_on1= 5.5e5;
%     k_off1 = 7.6e-3;
%     k_on2 = 7.57e5;
%     k_off2 = 3.3e-34;
%     scale = 1.37e8;
%     B0 = 2.5e-11;%smax
% 
%     A0 = 7e-9; %concentration
% 
%     maxTime = max(timePoints);
% 
%     t_off = 1980;
    %%
    timePoints = time;
    %%%This is the deltaT for trying to numerically integrate the equation...
    maxTime = max(timePoints);
    deltaT = 0.2;
    
    timeOffStep = round(t_off/deltaT);
    internalTime = 0:deltaT:maxTime;

    numSteps = round(maxTime/deltaT) ;  %The number of steps to take
    x = [0];
    y = [0];
    %%
    %x.append(0) %%%Define initial conditions == 0
    %y.append(0)


    %%The numerical integration step....
    for n = 2:numSteps+1
        temp1 = (2*k_on1*A0*(B0-x(n-1)-y(n-1)) - k_off1*x(n-1)-k_on2*x(n-1)*A0 + 2*k_off2*y(n-1))*deltaT + x(n-1);
        x = [x temp1];
        temp2 = (k_on2*x(n-1)*A0 - 2*k_off2*y(n-1))*deltaT + y(n-1);
        y = [y temp2];

        if n >= timeOffStep
            A0 = 0;
        end


    end
    %%
    %tempResult = [(x+y) for x,y in zip(x,y)]  %This is the total amount of increased thickness at the spot...Sums up the arrays component by component

    tempResult = x+y;

    %%Now, we try to find the data points from the numerical integration that are closest to the time points from the real dataset.
    newResult = [];
    %%
    
    for i = 1:length(timePoints)
        
        j = timePoints(i);
        try
            index = find(abs(j-internalTime)<0.01,1);
            newResult = [newResult scale*tempResult(index)];

        catch

        end
    end
    %%Defining the last data point as == (n-1)th data point

    newResult = [newResult newResult(end)];

end



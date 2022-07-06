

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%IDENTICAL SPOTS VERSION

%%%%%%%%%%%%%%%%%%%%%%%%%%

Data = IRIS_getData(); % <--- 
rawData = IRIS_getData(); % <--- 

pixel_area = (4.5e-6)^2; %m
pixel_area = pixel_area/4; %effective
fullWell = 10867;
totalBit = 65536;



clear background
clear differential
clear spots
clear total_pixel
%clear Data
clear n_spots
%Data = input('Name of Spots to Analyze: '); %<--
n_spots = input('Insert the number of Spot ROIs: '); % <-- # of Spot ROIs in array



spot = zeros(n_spots,size(Data,2));
background = spot;
differential = spot;

spotBrightness = spot;

for i = 1:n_spots
    
    background(i,:) = Data((6+8*(i-1)),:);
    spot(i,:) = Data((2+8*(i-1)),:);
    differential(i,:) = spot(i,:) - background(i,:);
    
    total_pixel(i) = rawData((1+8*(i-1)),1)./pixel_area;
    spotBrightness(i,:) = rawData((2+8*(i-1)),:);
end
time = (1:size(Data,2));
%%

%%
%differential = differential/1.132; %%convert to thickness nm/100
differential = differential*10;%%convert to pm
differential = differential*1.3;%%convert to pg/mm2
%%
differential = differential - differential(:,1);
figure
plot(1:size(differential,2),differential);
%%

%%
% time = (1:size(Data,2));
% differential = differential - (1:size(differential,2))*0.035;
% figure
% plot(1:size(differential,2),differential);
%%
%FIX OUTLIERS/EVENTS

%% regular method
%figure
%b = filloutliers(norm,'previous','percentiles',[1 90]);
norm = filloutliers(differential,'previous',2);
%plot(1:size(differential,2),differential)
%% mean method
%figure
%b = filloutliers(b,'previous','percentiles',[1 90]);
norm = filloutliers(norm,'previous','mean');
norm(isnan(norm)) = 0;
%plot(1:size(norm,2),norm)

%%

figure
meanc = sum(norm,1)/n_spots;
plot(1:size(norm,2),meanc)
%%

%Here you can take a flat portion or take the full data
%This will not make a difference since we are fitting the data not
%averaging for 'noiseless curve'

%Flat portion in electrons
flat = (input('Insert the start of Flat portion: '):input('Insert the end of Flat portion: '));
flat_portion = norm(:,flat); %take any interval you want

flat_brightness = spotBrightness(:,flat);
%%
%OPTIONAL: Convert it to electrons/area, you will later on multiply by area
% flat_portion = fullWell.*flat_portion./totalBit; %convert to electrons
% plot(1:length(flat_portion),flat_portion); % in electrons

%%

%%
%Here we analyze the SNR of identical spots, averaging them 1 to all
%Combinations of averaged spots from 2 to all-1 are recorded
%Conbination sizes exceed 10000, instead of choosing combinations I randomized the selections n times

%FOR LARGE NUMBER OF SPOTS
clear SNR
clear SNR_snr
clear total_e
clear total_e_foo
clear total_e_eachSpot

for i = 1:size(flat_portion,1)
    
    i
    if i == 1
        filtered = smoothdata(flat_portion,2,'sgolay',21); 
        diff = (flat_portion - filtered).^2;
        STD = sqrt(sum(diff,2)); STD = STD./sqrt(length(diff));
        
        %FOR differential SNR (this makes more sense in bio point of veiew)
        SNR{i} = median(filtered,2)./STD; 
        %SNR{i} = flat_portion(:,1)./STD; 
        
        %FOR spot brightness SNR (allison does it like this)
%         SNR{i} = mean(spot(:,flat),2)./STD;
         

        total_e_eachSpot = total_pixel.*mean(flat_brightness,2).'*fullWell/totalBit;
         
        
        total_e{i} = total_e_eachSpot;
        plot(1:size(filtered,2),filtered);hold on;plot(1:size(flat_portion,2),flat_portion);
        
    elseif i == size(flat_portion,1)
    
        averaged = sum(flat_portion,1)./n_spots;
        filtered = smoothdata(averaged,'sgolay',21); 
        diff = (averaged - filtered).^2;
        STD = sqrt(sum(diff,2)); STD = STD./sqrt(length(diff));
        
        %FOR differential SNR
        SNR{i} = mean(filtered,2)./STD;
        SNR{i} = filtered(1)./STD; 
        %FOR spot brightness SNR
%         SNR{i} = mean(sum(spot(:,flat),1)/n_spots,2)./STD;
        
        total_e{i} = sum(total_e_eachSpot(:));
        
    else
        for n = 1:99 %Change according to your ROI set size
            
            ind = randperm(size(flat_portion,1));
            randomMat = flat_portion(ind, :);
            randomActual = flat_brightness(ind, :);
            randomE = total_e_eachSpot(ind);
            
            averaged = sum(randomMat((1:i),:),1)./i; 
            actualLevel = sum(randomActual((1:i),:),1)./i; 
            
            filtered = smoothdata(averaged,2,'sgolay',21); 
            diff = (averaged - filtered).^2;
            STD = sqrt(sum(diff,2)); STD = STD./sqrt(length(diff));
            
            %FOR differential SNR
            SNR_snr(n,:) = mean(filtered,2)./STD;
            SNR_snr(n,:) = filtered(1)./STD; 
            %FOR spot brightness SNR
%             SNR_snr(n,:) = mean(actualLevel,2)./STD;
            
            total_e_foo(n) = sum(randomE(1:i));
        end
        SNR{i} = SNR_snr;
        total_e{i} = total_e_foo;
    end

    
end

%%
%%Get the median SNR
clear averaged_e
figure
medians = [];mads = [];madsX = [];

for i = 1:numel(SNR)
   
    medians = [medians median(SNR{i})];mads = [mads mad(SNR{i})];madsX = [madsX mad(total_e{i})];
    ind = find(SNR{i}==median(SNR{i}));
    averaged_e(i) = mean(total_e{i}(ind));
   
end
y = medians; err = mads; errX = madsX;
scatter(averaged_e,y)
ylabel('SNR');
xlabel('# of Electrons')%<---- is this right? 


%%
%Histograms
bin = 20; %change this for histograms
figure
if numel(SNR)>=326, k = 326; 
    else, k = numel(SNR);
end
    
for i = 1:k 
    
    subplot(7,3,i); %change subplot settings according to numel(SNR)
    histogram(SNR{i},bin) 
    title(['#', num2str(i), ' Median-MAD ', num2str(median(SNR{i})),'-',num2str(mad(SNR{i}))])
    
end

%%
%SNR = SNR(1:150)
y(1) = []; averaged_e(1) = []; err(1)=[];
%%
%averaged_e = total_pixel(1)*pixel_area*(1:length(y));
[fitresult, gof] = createFit1(averaged_e, y);

[fitresult2, gof2] = createFit3(averaged_e, y);

%%
axes = linspace(averaged_e(1),averaged_e(end),100);
coeffvals= coeffvalues(fitresult);coeffvals2= coeffvalues(fitresult2);
a1 = coeffvals(1);a2 = coeffvals2(1);
b1 = 0.5; b2 = coeffvals2(2);

figure
plot(axes,a1*(axes).^b1,'g','LineWidth',2);hold on; 
plot(axes,a2*(axes).^b2,'r','LineWidth',2);hold on;


%e = errorbar(averaged_e,y,err,'o');
e = errorbar(averaged_e(1:25:end),y(1:25:end),err(1:25:end),'o','MarkerSize',4,...
    'MarkerEdgeColor','blue','MarkerFaceColor','blue');
hold on;
e.Color = 'b';


xlabel( 'Number of Electrons' );ylabel( 'SNR' );


grid on
title('Identical Spots SNR Fit: a*(x)^b')
legend(['b = 0.5 R^2 = ' num2str(gof.rsquare)], ['b = ' num2str(b2) ' R^2 = ' num2str(gof2.rsquare)],'Data','Location','SouthEast');
%%

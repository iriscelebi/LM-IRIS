
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
%REMEMBER TO CHANGE INDEXES OF SPOTS MANUALLY
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%clc


%HANDSPOT VERSION

%%%%%%%%%%%%%%%%%%%%%%%%%%
thicknessData = IRIS_getData(); % <--- 

rawData = IRIS_getData(); % <--- 


pixel_area = (4.5e-6)^2; %m
pixel_area = pixel_area/4; %effective
fullWell = 10867;
totalBit = 65536;
%%
clear background
clear spots
clear total_pixel
%clear Data
clear n_spots
%Data = input('Name of Spots to Analyze: '); %<--
n_spots = input('Insert the number of Spot ROIs: '); % <-- # of Spot ROIs in array


spot = zeros(n_spots,size(thicknessData,2));
background = spot;
differential = spot;
total_pixel = spot;
spotBrightness = spot;

for i = 1:n_spots
    
    background(i,:) = thicknessData((6+8*(i-1)),:);
    spot(i,:) = thicknessData((2+8*(i-1)),:);
    differential(i,:) = spot(i,:) - background(i,:);
    
    total_pixel(i) = rawData((1+8*(i-1)),1)./pixel_area;
    spotBrightness(i,:) = rawData((2+8*(i-1)),:);
end
time = (1:size(thicknessData,2));
%%
%differential = differential/1.32; %%convert to thickness nm/100
differential = differential*10;%%convert to pm
differential = differential*1.3;%%convert to pg/mm2


%%
% figure
% plot(1:size(spot,2),differential);
% trendData = sum(differential(:,(270:350)),1)./26;
%%
% 
% slope = -0.187;
% time = 1:size(differential,2);
% differential = differential(:,:) - repmat(slope*time,26,1);
% figure
% plot(time,differential)



%%
%FIX OUTLIERS/EVENTS
norm = differential - differential(:,1);
figure
plot(1:size(differential,2),norm);
%% regular method
figure
%b = filloutliers(b,'previous','percentiles',[1 90]);
norm = filloutliers(norm,'previous',2);
norm = filloutliers(norm,'previous','mean');
plot(1:size(norm,2),norm)

%%

%Here you can take a flat portion or take the full data
%This will not make a difference since we are fitting the data not
%averaging for 'noiseless curve'


%Flat portion in electrons
flat = (input('Insert the start of Flat portion: '):input('Insert the end of Flat portion: '));
flat_portion = norm(:,flat); %take any interval you want

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%flat_portion = flat_portion./repmat(median(flat_portion,2),1,size(flat_portion,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flat_brightness = spotBrightness(:,flat);

% %%
% %OPTIONAL: Convert it to electrons/area, you will later on multiply by area
% flat_portion = fullWell.*flat_portion./totalBit; %convert to electrons
% plot(1:length(flat_portion),flat_portion); % in electrons


ind1 = 1:4;ind10 = 5:8;ind30 = [9 10 11 12];ind100 = [13 14 15 16];ind300 = [17 18 19 20];indH = [21 22 23];

spot300 = flat_portion(ind300,:);pixel300 = total_pixel(ind300);brightness300 = flat_brightness(ind300,:);
spot100 = flat_portion(ind100,:);pixel100 = total_pixel(ind100);brightness100 = flat_brightness(ind100,:);
spot30 = flat_portion(ind30,:);pixel30 = total_pixel(ind30);brightness30 = flat_brightness(ind30,:);
spot10 = flat_portion(ind10,:);pixel10 = total_pixel(ind10);brightness10 = flat_brightness(ind10,:);
spot1 = flat_portion(ind1,:);pixel1 = total_pixel(ind1);brightness1 = flat_brightness(ind1,:);
spothand = flat_portion(indH,:);pixelhand = total_pixel(indH);brightnessHand = flat_brightness(indH,:);



%%
spot30(4,:) = [];
spot100(3,:) = []; 
spot300(3,:) = [];
spothand(1,:) = [];
%%
%Average them 1 to all

%FOR HAND SPOTS

SNRall = []; total_e_all = [];
allSpots{1} = spot1; allSpots{2} = spot10; allSpots{3} = spot30; allSpots{4} = spot100; allSpots{5} = spot300; allSpots{6} = spothand; 
allPixels{1} = pixel1; allPixels{2} = pixel10; allPixels{3} = pixel30; allPixels{4} = pixel100; allPixels{5} = pixel300; allPixels{6} = pixelhand; 
allBrightness{1} = brightness1;allBrightness{2} = brightness10;allBrightness{1} = brightness1;allBrightness{3} = brightness30;allBrightness{4} = brightness100;allBrightness{5} = brightness300;allBrightness{6} = brightnessHand;

ind1 = 1:size(spot1,1);ind10 = ind1(end)+1:size(spot10,1)+ind1(end);ind30 = ind10(end)+1:size(spot30,1)+ind10(end);ind100 = ind30(end)+1:size(spot100,1)+ind30(end);ind300 = ind100(end)+1:size(spot300,1)+ind100(end);indH = ind300(end)+1:size(spothand,1)+ind300(end);


time = (1:size(spot1,2));

for j = 1:numel(allSpots)
   
    spot_type = allSpots{j};
    pixels = allPixels{j};
    brightness = allBrightness{j};
   
    for i = 1:size(spot_type,1)

        ind = nchoosek(1:size(spot_type,1),i);
        figure
        clear averaged    
        clear pixel_sum   
        clear averaged_brightness
       
        for n = 1:size(ind,1)
            n
            averaged(n,:) = sum(spot_type(ind(n,:),:),1)./size(ind,2);
            averaged_brightness(n,:) = sum(brightness(ind(n,:),:),1)./size(ind,2);
            pixel_sum(n) = sum(pixels(ind(n,:)));
        end
        
            filtered = smoothdata(averaged,2,'sgolay',21); 
            diff = (averaged - filtered).^2;
            STD = sqrt(sum(diff,2)); STD = STD./sqrt(length(diff));
            SNR{i} = mean(filtered,2)./STD;
            
            total_e{i} = pixel_sum.*mean(averaged_brightness,2).'*fullWell/totalBit; %*electron_eachPixel;
        

            plot(time,filtered); hold on;plot(time,averaged);hold on

    end
    SNRall = [SNRall SNR];
    total_e_all = [total_e_all total_e];
    clear SNR
    clear total_e
end



%%

%%TAKE OUT DATA IF YOU WANT

SNRall = [SNRall(1:10) SNRall(end-2:end)];
total_e_all = [total_e_all(1:10) total_e_all(end-2:end)];


%%
%%Get the median SNR

%figure
medians = [];mads = [];madsX = [];
clear averaged_e
for i = 1:numel(SNRall)
   
    medians = [medians median(SNRall{i})];mads = [mads mad(SNRall{i})];
    
    %ind = find(SNRall{i}==median(SNRall{i}));
    [foo,ind] = min(abs(SNRall{i}-median(SNRall{i})));
    averaged_e(i) = mean(total_e_all{i}(ind));
   
end
y = medians; err = mads; 
SNRmedian = y;

%%
%Histograms
% bin = 20; %change this for histograms
% figure
% for i = 1:numel(SNR)
%     
%     subplot(7,3,i); %change subplot settings according to numel(SNR)
%     histogram(SNR{i},bin) 
%     title(['#', num2str(i), ' Median-MAD ', num2str(median(SNR{i})),'-',num2str(mad(SNR{i}))])
%     
% end

%%
%SORT THE VALUES for HandSpot chip

%%AREA VERSION
%%Here you will need to enter the area for one spot of every type
%%% ex. area1 = area of one spot of one drop 
%%% area10 = area of one spot ten drop

axes1 = pixel1(1)*pixel_area*(1:size(spot1,1));
axes10 = pixel10(1)*pixel_area*(1:size(spot10,1));
axes30 = pixel30(1)*pixel_area*(1:size(spot30,1));
axes100 = pixel100(1)*pixel_area*(1:size(spot100,1));
axes300 = pixel300(1)*pixel_area*(1:size(spot300,1));
axesH = pixelhand(1)*pixel_area*(1:size(spothand,1));

areas = [axes1 axes10 axes30 axes100 axes300 axesH];

%%
[areas,I] = sort(areas);
Xaxes = areas;

%%This is to get rid of areas greater than 1mm
% % ind = find(areas>1.1);
% % areas(ind) = [];SNR_sort(ind) = [];


%%

%%TOTAL E VERSION

axes1 = averaged_e(ind1);
axes10 = averaged_e(ind10);
axes30 = averaged_e(ind30);
axes100 = averaged_e(ind100);
axes300 = averaged_e(ind300);
axesH = averaged_e(indH);
%%
[sorted_electron,I] = sort(averaged_e);
Xaxes = sorted_electron;



%%
%%Check if this part is correct


SNR1 = SNRmedian(ind1);
SNR10 = SNRmedian(ind10);
SNR30 = SNRmedian(ind30);
SNR100 = SNRmedian(ind100);
SNR300 = SNRmedian(ind300);
SNRH = SNRmedian(indH);
%%
SNRmedian = SNRmedian(I);
%%
%%FIT

[fitresult, gof] = createFit1(Xaxes, SNRmedian);

[fitresult2, gof2] = createFit3(Xaxes, SNRmedian);

%%
%Create Plots
figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%I seperated them to have different colors, they will have different names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear plotHandles 
clear plotLabels
plotHandles(1) = errorbar(axes1,SNR1,err(ind1),'o','MarkerFaceColor','k','Color','k','MarkerSize',8);hold on
plotHandles(2) = errorbar(axes10,SNR10,err(ind10),'o','MarkerFaceColor','r','Color','r','MarkerSize',8);hold on
plotHandles(3) = errorbar(axes30,SNR30,err(ind30),'o','MarkerFaceColor','g','Color','g','MarkerSize',8);hold on
plotHandles(4) = errorbar(axes100,SNR100,err(ind100),'o','MarkerFaceColor','b','Color','b','MarkerSize',8);hold on
plotHandles(5) = errorbar(axes300,SNR300,err(ind300),'o','MarkerFaceColor','m','Color','m','MarkerSize',8);hold on
plotHandles(6) = errorbar(axesH,SNRH,err(indH),'o','MarkerFaceColor','c','Color','c','MarkerSize',8);hold on
%%
axes_plot = linspace(Xaxes(1),Xaxes(end),1000);
coeffvals= coeffvalues(fitresult);coeffvals2= coeffvalues(fitresult2);
a1 = coeffvals(1);a2 = coeffvals2(1);
b1 = 0.5; b2 = coeffvals2(2);
p1 = plot(axes_plot,a1*(axes_plot).^b1,'Color',[0.8500 0.3250 0.0980],'LineWidth',1);hold on;
p2 = plot(axes_plot,a2*(axes_plot).^b2,'Color',[0.9290 0.6940 0.1250],'LineWidth',1);


plotLabels{1} = '1 Drop';
plotLabels{2} = '10 Drop';
plotLabels{3} = '30 Drop';
plotLabels{4} = '100 Drop';
plotLabels{5} = '300 Drop';
plotLabels{6} = 'HandSpot';

plotHandles(7) = p1;
plotLabels{7} = ['b = 0.5 R^2 = ' num2str(gof.rsquare)]; 
plotHandles(8) = p2;
plotLabels{8} = ['b = ' num2str(b2) ' R^2 = ' num2str(gof2.rsquare)];


grid on
title('SNR Fit: a*(x)^b')

legend(plotHandles, plotLabels,'Location','SouthEast');
ylabel('SNR');xlabel('Area m^2')
%ylabel('SNR');xlabel('Number of Electrons')


%%
figure
axes_plot = linspace(xfoo(1),xfoo(end),1000);
coeffvals= coeffvalues(fitresult);coeffvals2= coeffvalues(fitresult2);
a1 = coeffvals(1);a2 = coeffvals2(1);
b1 = 0.5; b2 = coeffvals2(2);
p1 = plot(axes_plot,a1*(axes_plot).^b1,'Color',[0.8500 0.3250 0.0980],'LineWidth',1);hold on;
p2 = plot(axes_plot,a2*(axes_plot).^b2,'Color',[0.9290 0.6940 0.1250],'LineWidth',1);

legend(['b = 0.5 R^2 = ' num2str(gof.rsquare)],['b = ' num2str(b2) ' R^2 = ' num2str(gof2.rsquare)])
%%
figure
plotHandles(1) = errorbar(axes1(1),SNR1(1),err(ind1(1)),'o','MarkerFaceColor','k','Color','k','MarkerSize',8);hold on
plotHandles(2) = errorbar(axes10(1),SNR10(1),err(ind10(1)),'o','MarkerFaceColor','r','Color','r','MarkerSize',8);hold on
plotHandles(3) = errorbar(axes30(1),SNR30(1),err(ind30(1)),'o','MarkerFaceColor','g','Color','g','MarkerSize',8);hold on
plotHandles(4) = errorbar(axes100(1),SNR100(1),err(ind100(1)),'o','MarkerFaceColor','b','Color','b','MarkerSize',8);hold on
plotHandles(5) = errorbar(axes300(1),SNR300(1),err(ind300(1)),'o','MarkerFaceColor','m','Color','m','MarkerSize',8);hold on
%plotHandles(6) = errorbar(axesH(1),SNRH(1),err(indH(1)),'o','MarkerFaceColor','c','Color','c','MarkerSize',8);hold on

%%
scatter([1537187702.99119 1078906657.92243],[216.430055200354 192.399593472462],60,'filled')
%%
foohandaxis = [1078906657.92243 1537187702.99119 ];foohandsnr = [192.399593472462 216.430055200354 ];
xfoo = [axes1(1) axes10(1) axes30(1) axes100(1) axes300(1) foohandaxis];
yfoo = [SNR1(1) SNR10(1) SNR30(1) SNR100(1) SNR300(1) foohandsnr];
[fitresult, gof] = createFit1(xfoo, yfoo);

[fitresult2, gof2] = createFit3(xfoo, yfoo);

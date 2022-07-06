
%% Fit exponentials
%SNR = SNR_hand;
xaxis = 1:numel(y);
[fitresult, gof] = createFit1(xaxis, y);
[fitresult2, gof2] = createFit2(xaxis, y);


%%
axes = linspace(1,numel(y),100);
coeffvals= coeffvalues(fitresult);coeffvals2= coeffvalues(fitresult2);
a1 = coeffvals(1);a2 = coeffvals2(1);
b1 = 0.5; b2 = coeffvals2(2);

figure

scatter(1:numel(y),y,80,'filled','b');hold on; 
plot(axes,a1*(axes).^b1,'g','LineWidth',2);hold on; 
plot(axes,a2*(axes).^b2,'r','LineWidth',2);
xlabel( 'number of spots' );ylabel( 'SNR' );
e = errorbar(1:numel(y),y,err,'o');hold on;
e.Color = 'b';
grid on
title('Epoxy SNR Fit: a*(x)^b')
legend('Data',['b = 0.5 R^2 = ' num2str(gof.rsquare)], ['b = ' num2str(b2) ' R^2 = ' num2str(gof2.rsquare)],'Location','SouthEast');


%%
figure
% 
% scatter(area1,SNR_one,100,'filled');hold on
% scatter(area2,SNR_ten,100,'filled');hold on
% scatter(area3,SNR_thirty,100,'filled');hold on
% scatter(area4,SNR_hundred,100,'filled');hold on
% scatter(area5,SNR_threeHundred,100,'filled');hold on
% scatter(areaH(1),SNR_handSpot,100,'filled');hold on

scatter(area1,SNR_oneDrop,100,'filled');hold on
scatter(area2,SNR_tenDrop,100,'filled');hold on
scatter(area3,SNR_thirtyDrop,100,'filled');hold on
scatter(area4,SNR_hundredDrop,100,'filled');hold on
scatter(area5,SNR_threeHundredDrop,100,'filled');hold on
scatter(areaH(1:3),SNR_handSpot,100,'filled');hold on

%plot(axes,a1*(axes).^b1+d1,'LineWidth',2);hold on; plot(axes,a2*(axes).^b2,'LineWidth',2);

xlabel( 'Area (m^2)' );ylabel( 'SNR' );
grid on
title('Epoxy: SNR vs Averaged Area')
%legend('Data',['b = 0.5 R^2 = ' num2str(gof.rsquare)], ['b = ' num2str(b2) ' R^2 = ' num2str(gof2.rsquare)], 'Location', 'SouthEast','FontSize',14);
legend('1 drop','10 drop','30 drop','100 drop','300 drop','Hand Spot')
%legend('1 drop','10 drop','30 drop','100 drop','300 drop','Hand Spot',['b = 0.5 R^2 = ' num2str(gof.rsquare)], ['b = ' num2str(b2) ' R^2 = ' num2str(gof2.rsquare)])
%%
areas = [area1 area2 area3 area4 area5 areaH(1)];
SNRs = [SNR_one SNR_ten SNR_thirty SNR_hundred SNR_threeHundred SNR_hand];
%SNRs = [SNR_oneDrop SNR_tenDrop SNR_thirtyDrop SNR_hundredDrop SNR_threeHundredDrop SNR_handSpot];

[areas,I] = sort(areas);
SNR_sort = SNRs(I);

[fitresult, gof] = createFit(areas, SNR_sort);
[fitresult2, gof2] = snr_powerfit2(areas, SNR_sort);

axes = linspace(areas(1),areas(end),1000);
coeffvals= coeffvalues(fitresult);coeffvals2= coeffvalues(fitresult2);
a1 = coeffvals(1);a2 = coeffvals2(1);
b1 = 0.5; b2 = coeffvals2(2);
d1 = coeffvals(2);
%c1 = coeffvals(2);c2 = coeffvals2(3);
figure
scatter(areas,SNRs,'filled');hold on; plot(axes,a1*(axes).^b1+d1,'LineWidth',2);hold on; plot(axes,a2*(axes).^b2,'LineWidth',2);
xlabel( 'number of spots' );ylabel( 'SNR' );
grid on
title('SNR Fit: a*(x)^b')
legend('Data',['b = 0.5 R^2 = ' num2str(gof.rsquare)], ['b = ' num2str(b2) ' R^2 = ' num2str(gof2.rsquare)], 'Location', 'SouthEast','FontSize',14);

%%
%For plotting all curves in one 
figure
STDarray = zeros(1,length(averages_forall(:,1)));
maxes = averages_forall(:,200);
normalized = averages_forall./(repmat(maxes,1,length(y(1,:))));

for i = 1: length(averages_forall(:,1))
    y = normalized(i,:);
    thick = i*0.5;
    plot(time,y,'LineWidth',thick);hold on
end
legend('1 drop','10 drop','30 drop','100 drop','300 drop','Hand Spot')



 
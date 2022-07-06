

%%



%%
%%%%%%%%
%%
%%
%SORT THE VALUES for Matts chip
%%
m = 1e6;
areas = [spot1; spot10; spot30; spot100; spot300; spothand; spotBSA];
areas = areas*m;
SNRmedian = median(areas);
[areas,I] = sort(areas);SNR_sort = SNRmedian(I); %sorts the arrayed areas from min to max (sorting area by size of area, combined and non-combined)

med = median(I);SNR_sort = areas(med); % (randi([1,1])
%areaSmall = [area1 area2 area3 area4 areaH];areaSmall = areaSmall*m;
%[areaSmall,I] = sort(areaSmall);SNR_sort = SNRmedian(I);

%%
ind = find(areas>1.1);
%areas(ind) = [];SNR_sort(ind) = [];


%%
[fitresult, gof] = createFit1(areas, SNR_sort);
[fitresult2, gof2] = createFit2(areas, SNR_sort);

%%
axes = linspace(areas(1),areas(end),1000);
coeffvals= coeffvalues(fitresult);coeffvals2= coeffvalues(fitresult2);
a1 = coeffvals(1);a2 = coeffvals2(1);
b1 = 0.5; b2 = coeffvals2(2);


p1 = plot(axes,a1*(axes).^b1,'Color',[0.8500 0.3250 0.0980],'LineWidth',1);hold on;
p2 = plot(axes,a2*(axes).^b2,'Color',[0.9290 0.6940 0.1250],'LineWidth',1);
xlabel( 'Area' );ylabel( 'SNR' );


plotHandles(6) = p1;plotHandles(7) = p2;
plotLabels{6} = ['b = 0.5 R^2 = ' num2str(gof.rsquare)]; 
plotLabels{7} = ['b = ' num2str(b2) ' R^2 = ' num2str(gof2.rsquare)];


grid on
title('SNR Fit: a*(x)^b')
%legend('Data',['b = 0.5 R^2 = ' num2str(gof.rsquare)], ['b = ' num2str(b2) ' R^2 = ' num2str(gof2.rsquare)], 'Location', 'SouthEast','FontSize',14);
legend(plotHandles, plotLabels,'Location','SouthEast');
ylabel('SNR');xlabel('Area (mm^2)')

%%
%%
%%%%%%%%
%%
%%
%FOR Elis chip
%Histograms
means = [];stds = [];medians = [];mads = [];
figure
for i = 1:numel(SNR)
    
    subplot(7,3,i);histogram(SNR{i},20)
    means = [means mean(SNR{i})];stds = [stds std(SNR{i})];
    medians = [medians median(SNR{i})];mads = [mads mad(SNR{i})];
    title(['#', num2str(i), ' Median-MAD ', num2str(median(SNR{i})),'-',num2str(mad(SNR{i}))])
    
end
y = medians; err = mads;
%%
figure
plot(1:numel(means),means,1:numel(medians),medians);
figure
plot(1:numel(stds),stds,1:numel(mads),mads);
%%

figure
%errorbar(x,y,yneg,ypos,'o')
errorbar(1:numel(y),y,err,'o')
title('18 Spot SNR fit: a*x^b')
%%
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
title('18 Spot SNR Fit: a*(x)^b')
legend('Data',['b = 0.5 R^2 = ' num2str(gof.rsquare)], ['b = ' num2str(b2) ' R^2 = ' num2str(gof2.rsquare)],'Location','SouthEast');
%%


figure

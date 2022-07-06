% %%MOCK DATA
% k_on = 2.28e4; 
% k_off = 2.36e-05;
% 
% smax =  1;
% scale =  100;
% 
% t = 1:6:2500;
% stop_time = 1309;
% stop_time_ind = find(t == stop_time);
% [y, assoc, dissoc] = langmuirModel(k_on, k_off, smax, scale, t, stop_time);
% t_assoc = t(1:stop_time_ind);

%REAL DATA
injection = 13;
assoc = signal_LAC1_10ugmL(injection:end);
dissoc = signal_LAC1_PBS_10ugmL;
y = [assoc dissoc];

t_assoc = 1:6:6*size(assoc,2); t = 1:6:6*size(y,2);

assoc = assoc/max(y);dissoc = dissoc/max(y); %normalize
y = [assoc dissoc];
scatter(1:size(y,2),y)
%%
%%FIT TO LANGMUIR
concentration = 70e-9; scale = 2;
%[fitresult, gof] = FullLangmuirFit(t_assoc, assoc, concentration, scale);
[fitresult, gof] = FullLangmuirFit(t, y, 70e-9, scale, 1000, 1300);

fitresult
coeffvals= coeffvalues(fitresult); 
a = coeffvals(1);koff = coeffvals(3);kon = coeffvals(4);

%%

y_est = langmuirModel(kon, koff, 1, a, t, t_assoc(end));
Rsq1 = 1 - sum((y - y_est).^2)/sum((y - mean(y)).^2);
figure; scatter(t, y);hold on;plot(t, y_est,'LineWidth',3)
legend(['Langmuir Fit R^2 = ' num2str(gof.rsquare)],'Data')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
% lambdaArray = 100:-5:5;
% 
% for i = 1:length(lambdaArray)
% 
%     lambda = lambdaArray(i);
%     noise =  poissrnd(lambda,size(t));
%     noisyData = noise+mocky; noisyData = noisyData - noisyData(1);
%     %figure;plot(noise);title('noise')
%     meanNoise(i) = sum(noise);
% 
% 
%     noisyData = noisyData - noisyData(1);
%     filtered = smoothdata(noisyData,'sgolay',41); 
%     figure;plot(t,noisyData);hold on;plot(t,filtered)
% 
%     diff = (noisyData - filtered).^2;
%     %figure;plot(sqrt(diff));title('diff')
%     STD = sqrt(sum(diff,2)); STD = STD./sqrt(length(diff));
%     SNR(i) = max(filtered)./STD; 
% end
% %%
% figure
% scatter(lambdaArray,SNR)
% 
% %%
% figure;
% scatter(lambdaArray,meanNoise)


%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number_of_curves = 20;
% noise =  poissrnd(50,size(time,2),number_of_curves);
% noisyData = repmat(mocky,number_of_curves,1)+noise.';
% noisyData = noisyData - noisyData(:,1);
% 
% axis = 1:number_of_curves;
% %%
% scale = 1.5;
% clear SNR
% for i = 1:number_of_curves
%     i
%     summed = sum(noisyData(1:i,:),1)/i;
%     filtered = smoothdata(summed,2,'sgolay',41); 
%     diff = (summed - filtered).^2;
%     STD = sqrt(sum(diff,2)); STD = STD./sqrt(length(diff));
%     SNR(i) = max(filtered)./STD;
%     
%     %figure; plot(t,summed);hold on;plot(t,filtered)
%     normSummed = summed/max(filtered);
%     normFiltered = filtered/max(filtered);
%     %[fitresult, gof] = LangmuirFit(t_assoc, normFiltered(1:188), concentration, scale);
%     [fitresult, gof] = FullLangmuirFit(time, normFiltered, 70e-9, scale, 1000, 1200);
% 
%     coeffvals= coeffvalues(fitresult); 
%     a(i) = coeffvals(1);koff(i) = coeffvals(3);kon(i) = coeffvals(4);
%     y_est = langmuirModel(kon(i), koff(i), 70e-9, a(i), time, stopTime);
%     Rsq1 = 1 - sum((normSummed - y_est).^2)/sum((normSummed - mean(normSummed)).^2);
%     
%     
%     figure; scatter(time, normSummed);hold on;plot(time, y_est,'LineWidth',4)
%     plot(time,normFiltered,'LineWidth',2)
%     legend(['Langmuir Fit R^2 = ' num2str(Rsq1)],'Data','Location','SouthEast')
%     
% end
% %%
%%
for i =1:6
    figure
    myData = firsts(i,70:360);
    time = 1:6:6*length(myData);
    filtered = smoothdata(myData,'sgolay',21);
    [fitresult, gof] = FullLangmuirFit(time, myData/max(myData), 70e-9, 1.5, 1000, 1200);
    coeffvals= coeffvalues(fitresult); 
    a = coeffvals(1);koff = coeffvals(3);kon = coeffvals(4);stopTime = coeffvals(5);

    y_est = langmuirModel(kon, koff, 7e-8, a, time, stopTime);
    plot(time, myData);hold on;plot(time, y_est*max(myData),'LineWidth',3);hold on;plot(time, filtered,'LineWidth',2)
    legend(['Langmuir Fit R^2 = ' num2str(gof.rsquare)],'Data')
    fitresult
end
figure
filtered = smoothdata(firsts(:,70:end),2,'sgolay',41);
plot(1:6:6*length(filtered(:,70:end)),filtered(:,70:end))
figure 
plot(1:6:6*length(firsts(:,70:end)),firsts(:,70:end))
legend('Single Droplet','10 Droplet','30 Droplet','100 Spot','300 Droplet','Hand Spot')

% mocky = langmuirModel(kon, koff, 7e-8, 100, time, stopTime);
% figure;
% plot(time, mocky,'LineWidth',3);hold on
%%
figure
%myData = fooData(:,70:360)/2;
time = 1:6:6*length(myData);
%filtered = smoothdata(fooData/2,'sgolay',11);
[fitresult, gof] = FullLangmuirFit(time, myData/max(myData), 70e-9, 1.5, 1000, 1200);
coeffvals= coeffvalues(fitresult); 
a = coeffvals(1);koff = coeffvals(3);kon = coeffvals(4);stopTime = coeffvals(5);

y_est = langmuirModel(kon, koff, 7e-8, a, time, stopTime);
scatter(1:6:6*length(myData), myData,80);hold on;plot(time, y_est*max(myData),'LineWidth',4);hold on;plot(1:6:6*length(filtered(70:end)), filtered(70:end),'LineWidth',2)
legend(['Langmuir Fit R^2 = ' num2str(gof.rsquare)],'Data','SG Filtered')
fitresult

%%
kon = 2.28e4; 
koff = 2.36e-05;
t = 1:6:2500;
stop_time = 1309;
mocky = langmuirModel(kon, koff, 7e-9, 1000, t, 1280);

figure;
plot(t, mocky,'LineWidth',3);hold on
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
number_of_curves = 201;
noise =  poissrnd(50,size(t,2),number_of_curves);
noisyData = repmat(mocky,number_of_curves,1)+noise.';
noisyData = noisyData - noisyData(:,1);

axis = 1:number_of_curves;
%%
%%%%%%%%%%%%%%%%%%%%%%
%FOR LARGE NUMBER OF SPOTS
clear SNR
clear SNR_snr
clear total_e
clear total_e_foo
clear total_e_eachSpot
clear a A
clear koff KOFF
clear kon KON
    

concentration = 7e-9;
scale = 2;
stopLower = 1150;
stopUpper = 1300;


for i = 1:10:number_of_curves
    i
    
    if i == 1
        filtered = smoothdata(noisyData,2,'sgolay',41); 
        diff = (noisyData - filtered).^2;
        STD = sqrt(sum(diff,2)); STD = STD./sqrt(length(diff));
        SNR{i} = max( filtered ,[], 2 )./STD;
        figure
        plot(1:size(filtered,2),filtered);hold on;plot(1:size(noisyData,2),noisyData);
        normAveraged = noisyData./repmat(max( noisyData ,[], 2 ),1,size(noisyData,2));
        
        
        for n = 1:number_of_curves
        n
        [fitresult, gof] = FullLangmuirFit(time, normAveraged(n,:), concentration, scale, stopLower, stopUpper);

        coeffvals= coeffvalues(fitresult); 
        a(n) = coeffvals(1);koff(n) = coeffvals(3);kon(n) = coeffvals(4);
        
        end
        A{i} = a;
        KOFF{i} = koff;KON{i} = kon;
        
    elseif i == size(noisyData,1)
    
        averaged = sum(noisyData,1)./number_of_curves;
        filtered = smoothdata(averaged,'sgolay',41); 
        diff = (averaged - filtered).^2;
        STD = sqrt(sum(diff,2)); STD = STD./sqrt(length(diff));
        
        %FOR differential SNR
        SNR{i} = max( filtered ,[], 2 )./STD;
        
        normAveraged = averaged/max(averaged);
        [fitresult, gof] = FullLangmuirFit(time, normAveraged, concentration, scale, stopLower, stopUpper);

        coeffvals= coeffvalues(fitresult); 
        A{i} = coeffvals(1);KOFF{i} = coeffvals(3);KON{i} = coeffvals(4);
        


        
    else
        for n = 1:99 %Change according to your ROI set size
            
            ind = randperm(size(noisyData,1));
            randomMat = noisyData(ind, :);
            
            averaged = sum(randomMat((1:i),:),1)./i; 
          
            filtered = smoothdata(averaged,2,'sgolay',41); 
            diff = (averaged - filtered).^2;
            STD = sqrt(sum(diff,2)); STD = STD./sqrt(length(diff));
            
            %FOR differential SNR
            SNR_snr(n,:) = max( filtered ,[], 2 )./STD;
            
            
            normAveraged = averaged/max(averaged);
            [fitresult, gof] = FullLangmuirFit(time, normAveraged, concentration, scale, stopLower, stopUpper);

            coeffvals= coeffvalues(fitresult); 
            a(n) = coeffvals(1);koff(n) = coeffvals(3);kon(n) = coeffvals(4);
            
    
            
        end
        SNR{i} = SNR_snr;
        
        A{i} = a;
        KOFF{i} = koff;KON{i} = kon;
        
        clear a
        clear koff
        clear kon
        clear SNR_snr
        
    end

    
end
%%
%%

%%
%%Get the median SNR
SNR = SNR(~cellfun(@isempty, SNR));
KON = KON(~cellfun(@isempty, KON));
KOFF = KOFF(~cellfun(@isempty, KOFF));
A = A(~cellfun(@isempty, A));
figure
medians = [];mads = [];madsX = [];perErrKON = []; perErrKOFF = [];

for i = 1:numel(SNR)
   
    medians = [medians median(SNR{i})];mads = [mads mad(SNR{i})];
    ind = find(SNR{i}==median(SNR{i}));
    
    medianA(i) = A{i}(ind);medianKOFF(i) = KOFF{i}(ind);medianKON(i) = KON{i}(ind);
    
    perErrKON{i} = 100*(KON{i}-kon)/kon; perErrKOFF{i} = 100*(KOFF{i}-koff)/koff; perErrKD{i} = 100*(KOFF{i}./KON{i}-kd)/kd;
    perErrKON_median(i) = abs(100*(medianKON(i)-kon)/kon);
    perErrKOFF_median(i) = abs(100*(medianKOFF(i)-koff)/koff);
    
    
    maxErrKON(i) = max(abs(perErrKON{i}));minErrKON(i) = min(abs(perErrKON{i}));
    maxErrKOFF(i) = max(abs(perErrKOFF{i}));minErrKOFF(i) = min(abs(perErrKOFF{i}));
    maxErrKD(i) = max(abs(perErrKD{i}));
    
    
end
y = medians; err = mads; errX = madsX;
scatter(1:numel(y),y)
ylabel('SNR');
%xlabel('curves')

%%

Xaxis =81:2:number_of_curves;
[fitresult, gof] = createFit1(Xaxis, y);
[fitresult2, gof2] = createFit3(Xaxis, y);

axes = linspace(Xaxis(1),Xaxis(end),100);
coeffvals= coeffvalues(fitresult);coeffvals2= coeffvalues(fitresult2);
a1 = coeffvals(1);a2 = coeffvals2(1);
b1 = 0.5; b2 = coeffvals2(2);

figure
plot(Xaxis,a1*(Xaxis).^b1,'g','LineWidth',2);hold on; 
plot(Xaxis,a2*(Xaxis).^b2,'r','LineWidth',2);hold on;
scatter(Xaxis,y,'filled')
legend(['b = 0.5 R^2 = ' num2str(gof.rsquare)], ['b = ' num2str(b2) ' R^2 = ' num2str(gof2.rsquare)],'Data','Location','SouthEast');

%%
bin = 20; %change this for histograms
figure

    
for i = 1:2:numel(KON)-11
    
    %subplot(7,3,i); %change subplot settings according to numel(SNR)
    histfit(KON{i},10); hold on
    %bar3(KON{i}); hold on
    %title(['#', num2str(i), ' Median-MAD ', num2str(median(SNR{i})),'-',num2str(mad(SNR{i}))])
    
end


%%


%%
[N,edges] = histcounts(KON{1},20);
clear N
centers = edges + (edges(2)-edges(1));centers = centers(1:end-1);
plotX = linspace(edges(1),centers(end),100);
%plotX = linspace(0.5e4,3e4,100);
figure
for ind = 1:2:numel(KON)-1
    i = ind*2;
    N(i,:) = (histcounts(KON{ind},edges)).';
    
    pd = fitdist(KON{ind}.','Normal');
    distrubution(i,:) = normpdf(plotX,pd.mu,pd.sigma)*pd.sigma*sqrt(2*pi)*max(N(i,:))*1.0;
%     dist = normpdf(plotX,pd.mu,pd.sigma);
%     distrubution(i,:) = dist/sum(dist);
    plot3(ones(1,100)*i,plotX,distrubution(i,:),'LineWidth',4);hold on
    clear dist
end
%%
figure
b = bar3(centers(1:end)+(edges(2)-edges(1)),N.');
set(b,'FaceAlpha',0.30)
remove_empty_bars(b)    

    

%%
fooB = bivalentFunc(10340,0.0002475,106479,0.00052023,25.09, 7e-8, 2.103, 1100,myTime);
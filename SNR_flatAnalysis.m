
%Data = IRIS_getData(); 


%%
pixel_area = 1.725^2; %um
fullWell = 6168;
totalBit = 65536;
%%
clear background
clear spots
clear total_pixel
n_spots = input('Insert the number of Spot ROIs:'); % <-- REPLACE ME: changes with each data set

%spot = zeros(n_spots,length(Data));
%background = zeros(n_spots,length(Data));
%differential = zeros(n_spots,length(Data));

spot = zeros(n_spots,size(Data,2));
background = zeros(n_spots,size(Data,2));
differential = zeros(n_spots,size(Data,2));

for i = 1:n_spots
    
    background(i,:) = Data((6+8*(i-1)),:);
    spot(i,:) = Data((2+8*(i-1)),:);
    differential(i,:) = spot(i,:) - background(i,:);
    total_pixel(i) = Data((1+8*(i-1)),1)./pixel_area;
end




%%

%%

%Flat portion in electrons
flat = (input('Insert the start of Flat portion:'):input('Insert the end of Flat portion:'));
flat_portion = differential(:,flat); %take any interval you want
flat_portion = fullWell.*flat_portion./totalBit; %convert to electrons
plot(1:length(flat_portion),flat_portion); % in electrons
%%

%FOR LARGE NUMBER OF SPOTS
clear SNR
clear SNR_snr
clear total_e
clear total_e_foo
clear total_e_eachSpot
for i = 1:size(flat_portion,1)
    
    i
    if i == 1
        filtered = smoothdata(flat_portion,2,'sgolay',41); 
        diff = (flat_portion - filtered).^2;
        STD = sqrt(sum(diff,2)); STD = STD./sqrt(length(diff));
        
        %SNR{i} = mean(filtered,2)./STD;
        SNR{i} = mean(spot(:,flat),2)./STD;
        
        %total_e_eachSpot = total_pixel.*mean(filtered,2).'; %*electron_eachPixel;
        total_e_eachSpot = total_pixel.*mean(spot(:,flat),2).';
        
        total_e{i} = total_e_eachSpot;
        plot(1:30,filtered);hold on;plot(1:30,flat_portion);
        
    elseif i == size(flat_portion,1)
    
        averaged = sum(flat_portion,1)./n_spots;
        filtered = smoothdata(averaged,'sgolay',41); 
        diff = (averaged - filtered).^2;
        STD = sqrt(sum(diff,2)); STD = STD./sqrt(length(diff));
        
        %SNR{i} = mean(filtered,2)./STD;
        SNR{i} = mean(sum(spot(:,flat),1)/n_spots,2)./STD;
        total_e{i} = sum(total_e_eachSpot(:));
    else
        for n = 1:9999
            
            ind = randperm(size(flat_portion,1));
            randomMat = flat_portion(ind, :);
            randomActual = spot(:,flat);randomActual = randomActual(ind, :);
            randomE = total_e_eachSpot(ind);
            
            averaged = sum(randomMat((1:i),:),1)./i; 
            actualLevel = sum(randomActual((1:i),:),1)./i; 
            
            filtered = smoothdata(averaged,2,'sgolay',41); 
            diff = (averaged - filtered).^2;
            STD = sqrt(sum(diff,2)); STD = STD./sqrt(length(diff));
            
            %SNR_snr(n,:) = mean(filtered,2)./STD;
            SNR_snr(n,:) = mean(actualLevel,2)./STD;
            total_e_foo(n) = sum(randomE(1:i));
        end
        SNR{i} = SNR_snr;
        total_e{i} = total_e_foo;
    end

    
end

%%

%%

%%Electron Only medians

figure
medians = [];mads = [];madsX = [];

for i = 1:numel(SNR)
   
    medians = [medians median(SNR{i})];mads = [mads mad(SNR{i})];madsX = [madsX mad(total_e{i})];
    ind = find(SNR{i}==median(SNR{i}));
    averaged_e(i) = mean(total_e{i}(ind));
   
end
y = medians; err = mads; errX = madsX;
scatter(averaged_e,y)

%%

[fitresult, gof] = createFit1(averaged_e, y);
[fitresult2, gof2] = createFit2(averaged_e, y);

%%
axes = linspace(averaged_e(1),averaged_e(end),100);
coeffvals= coeffvalues(fitresult);coeffvals2= coeffvalues(fitresult2);
a1 = coeffvals(1);a2 = coeffvals2(1);
b1 = 0.5; b2 = coeffvals2(2);

figure

plot(axes,a1*(axes).^b1,'g','LineWidth',2);hold on; 
plot(axes,a2*(axes).^b2,'r','LineWidth',2);
xlabel( 'number of total electrons' );ylabel( 'SNR' );

%e = errorbar(averaged_e,y,err,'o');
e = errorbar(averaged_e,y,err,'o','MarkerSize',8,...
    'MarkerEdgeColor','blue','MarkerFaceColor','blue');
hold on;
e.Color = 'b';

grid on
title('22 Spot SNR Fit: a*(x)^b')
legend(['b = 0.5 R^2 = ' num2str(gof.rsquare)], ['b = ' num2str(b2) ' R^2 = ' num2str(gof2.rsquare)],'Data','Location','SouthEast');
%%

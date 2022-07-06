
Data = IRIS_getData(); 

%%
n_spots = input('Insert the number of Spot ROIs:');
spot = zeros(n_spots,size(Data,2));
background = zeros(n_spots,size(Data,2));
differential = zeros(n_spots,size(Data,2));
for i = 1:n_spots
    
    background(i,:) = Data((6+8*(i-1)),:);
    spot(i,:) = Data((2+8*(i-1)),:);
    differential(i,:) = spot(i,:) - background(i,:);
    
end

%%
%plot(1:length(Data),spot)
plot(1:size(Data,2),spot)
grid on
%%
 
 for i = 1:n_spots
     foe(i,:) = Data((6+8*(i-1)),:);
     
 end
%%
 thickness = foo-foe;
%%

thickness(2,:) = [];

%%
initialThickness = thickness(:,1);
highestSignal = thickness(:,end); deltaThickness = highestSignal - initialThickness;
initial_delta = [initialThickness deltaThickness];
figure
scatter(initial_delta(:,1),initial_delta(:,2));
xlabel('Initial thickness');ylabel('Delta thickness');

%%
figure
plot(1:length(thickness),thickness)
%%
%correct the injection start 
injectionStart = 1 ;
thickness_start = thickness(:,(injectionStart:end));
time = 1:length(thickness)-injectionStart+1;
figure
plot(time,thickness_start)
%%
%START HERE
injectionStop = 220;

normalized = thickness_start - thickness_start(:,1);
%normalized = normalized(:,(1:270));

time = 1:length(normalized);
maxes = normalized(:,injectionStop);
normalized = normalized./(repmat(maxes,1,length(normalized(1,:))));
%%
figure
plot(time,normalized)
%%

%%

clear SNR
for i = 1:size(thickness,1)
    
    ind = nchoosek(1:size(thickness,1),i);
    figure
    for n = 1:size(ind,1)
        n
        
        averaged(n,:) = sum(normalized(ind(n,:),:),1)./size(ind,2);
        filtered = smoothdata(averaged,2,'sgolay',41); 
        diff = (averaged - filtered).^2;
        STD = sqrt(sum(diff,2)); STD = STD./sqrt(length(diff));
        SNR{i} = 1./STD;
        
        plot(time,filtered); hold on;plot(time,averaged);hold on

    
    end
    clear averaged
    
end
%%
clear SNR
for i = 1:size(normalized,1)
    
    
    if i == 1
        filtered = smoothdata(normalized,2,'sgolay',41); 
        diff = (normalized - filtered).^2;
        STD = sqrt(sum(diff,2)); STD = STD./sqrt(length(diff));
        SNR{i} = 1./STD;
        
    elseif i == size(normalized,1)
    
        averaged = sum(normalized,1)./18;
        filtered = smoothdata(averaged,'sgolay',41); 
        diff = (averaged - filtered).^2;
        STD = sqrt(sum(diff,2)); STD = STD./sqrt(length(diff));
        SNR{i} = 1./STD;
    else
        for n = 1:10000
            n
            ind = randperm(size(normalized,1));
            randomMat = normalized(ind, :);
            averaged = sum(randomMat((1:i),:),1)./i; 
            filtered = smoothdata(averaged,2,'sgolay',41); 
            diff = (averaged - filtered).^2;
            STD = sqrt(sum(diff,2)); STD = STD./sqrt(length(diff));
            SNR_snr(n,:) = 1./STD;

        end
        SNR{i} = SNR_snr;
    end

    
end
%%
SNR_hand = SNR;
%%

clear SNR
for i = 1:size(normalized,1)
   
    ind = nchoosek(1:size(thickness,1),i);
    
    for n = 1:size(ind,1)
        n
        
        averaged(n,:) = sum(normalized(ind(n,:),:),1)./size(ind,2);
        filtered = smoothdata(averaged,2,'sgolay',41); 
        diff = (averaged - filtered).^2;
        STD = sqrt(sum(diff,2)); STD = STD./sqrt(length(diff));
        SNR{i} = 1./STD;
     
    end
    clear averaged
    
end


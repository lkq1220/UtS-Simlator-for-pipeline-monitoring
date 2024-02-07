%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
clc;close all;clear; % reset all
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
tic
%% Gains and Pt are converted into linear form
Gr=(10.^((8)/10));      %GW 6.5 dBi antenna gain
%https://lora.readthedocs.io/en/latest/
Gt=(10.^((2.15)/10));      %2.15 dBi 

Freq_Band = 470e6;         % 470 MHz (frequency band China)
Band_Width=125e3;          %Band Width 
SpeedLight  = 3e8;     

Channels = 48;           
Pt = 10^(17/10)/1000;    
MonteCarlo = 1e2;          % 1e5 for the results in paper 

GW_num=1000;

[Grid, Idx, Center_found] = latlon2UTM('bus.csv',GW_num);  %Coordinates of Pipeline

Pipeline = length(Grid);
Sensors_per_Pipeline = 1;

figure
gscatter(Grid(:,1),Grid(:,2),Idx);
hold on
h(1)=plot(Center_found(:,1),Center_found(:,2),'kx','MarkerSize',14, 'LineWidth',2);

R_Time = 900;            
Period = 60.*60.*1000;     % 1 hour in milliseconds
D_SNR = 10.^([-17.5]./10);
ToA=[823.296];
Dutycycle=ToA./(1000*R_Time);
gamma_db = 6; % dB - SIR Threshold
gamma = 10^(gamma_db/10); % linear - SIR Threshold

%paramters for soil characteristics 
Clay_input=3.7;
VWC_input=0.132;
Depth=0.6;
[RealSoilDielectric, ImagSoilDielectric] = clc_die(Clay_input, VWC_input, Freq_Band);

Noise_Figure=6;   
NoisePower=(10.^((-174+Noise_Figure+10*log10(Band_Width))/10))/1000;

for i=1:length(D_SNR)
    for n=1:GW_num
        [i n]
        Node_distance_Idex= sort(sqrt((Grid(Idx==n,1)-Center_found(n,1)).^2+(Grid(Idx==n,2)-Center_found(n,2)).^2));
        for c=1:length(Node_distance_Idex)
            pr(n,c,:) = City_U2Aloss(RealSoilDielectric, ImagSoilDielectric,Depth,Node_distance_Idex(c),Freq_Band,MonteCarlo,Pt,Gt,Gr);
            P_SNR(i,n,c) = (sum((reshape(pr(n, c, :),1,[])./NoisePower) >= D_SNR(i)))/MonteCarlo;
            Distance(i,n,c)=Node_distance_Idex(c);
        end
    end
end

figure
for count=1:length(D_SNR)
    for i=1:GW_num
        
        D_tmp=reshape(Distance(count, i, :),1,[]);
        PSNR_tmp=reshape(P_SNR(count,i,:),1,[]);

        if count==1
            PSNRSF11_average(i) = mean(PSNR_tmp(1:length(D_tmp(D_tmp~=0))));
            n1 = plot(D_tmp(D_tmp~=0)/1000,PSNR_tmp(1:length(D_tmp(D_tmp~=0))),'r*');
        end
        hold on 
    end
end
grid on
ylabel('Success probability', 'Interpreter', 'Latex');
xlabel('Distance (km)', 'Interpreter', 'Latex');
legend([n1],{'$P_{SNR}$ (SF11)'},'Interpreter', 'Latex');
z = find(~isnan(PSNRSF11_average));
PSNRSF11_system = mean(PSNRSF11_average(z));

%% PSIR: Considering interfering signals 
for count=1:1
    for n=1:1:GW_num
        [count n]
        Node_distance_Idex= sort(sqrt((Grid(Idx==n,1)-Center_found(n,1)).^2+(Grid(Idx==n,2)-Center_found(n,2)).^2));
        for c=1:1:length(Node_distance_Idex)
            discarded = 0;     % To calculate destructive collision at a given distance(as function of elevation angle)
            Transmissions = (Period/(R_Time*1000))*length(Node_distance_Idex)*Sensors_per_Pipeline; 
            for m=1:MonteCarlo
                simultaneous=0;
        %% calculating the number of interfering signals
                Transmission_per_channel=round(Transmissions/Channels);  % Transmission per channel
                if Transmission_per_channel==0
                    Transmission_per_channel=1;
                end
                clear TX_stamp
                TX_stamp  = rand(1,round(Transmission_per_channel))*3600;   %Generating transmissions
               %TX_stamp  = rand(1,round(Total_Wind_Turbines/8)*6)*3600;
        
                Tsort=sort(TX_stamp);
                TimeStamp=Tsort;
                TimeStampEnd = TimeStamp + (ToA(count)/1000);       %End time of transmission of devices
                Transmit = randi(length(TimeStamp));                % Selecting one random device
                % Ts: transmission started
                Ts= TimeStamp(Transmit);
                % Tend: transmission ended
                Tend= Ts + (ToA(count)/1000);
            
                % Finding simultaneous transmission (collisions) in that time window 
                
                Starting = find(TimeStamp >= Ts & TimeStamp <= Tend);
            
                %Starting: Previous transmissions in that time window?
                Ending = find(TimeStampEnd >= Ts & TimeStampEnd <= Tend);
                %unique : find nodes other than the desired one in that window 
                simultaneous = length(unique([Starting Ending]))-1;          % colliding signals

                if simultaneous > 0                                  % if there is any simultaneous transmission (if collision)
              %% Generating location of interfering signals
			  
                    index = randi(length(Node_distance_Idex),simultaneous,1);
                    dPropogation=zeros(1,simultaneous); 
                    dPropogation=Node_distance_Idex(index);
    %               Lo = eta * 10*log10(4*pi*do*(1/wavelength));
                    clear L_int
                    clear L_Lin_int
                    clear pr_h_g_I
                    L_int=zeros(1,simultaneous);
                    L_Lin_int=zeros(1,simultaneous);
            
                    for track=1:length(dPropogation)
                        pr_h_g_I_tmp(track) = mean(City_U2Aloss(RealSoilDielectric, ImagSoilDielectric,Depth, dPropogation(track),Freq_Band,MonteCarlo,Pt,Gt,Gr));
                    end
                %% total interferance = sum of all the interfering signals 
                    pr_h_g_I = sum(pr_h_g_I_tmp);            

                %% Received power of desired signal
                   pr_h_g_D = mean(City_U2Aloss(RealSoilDielectric, ImagSoilDielectric,Depth, Node_distance_Idex(c),Freq_Band,MonteCarlo,Pt,Gt,Gr));

            %% Data Collisions and Spreading Factor Orthogonality
            
             %"Two packets with the same spreading factor arriving at the same time... 
             ...on the same channel might result in a collision. 
              ...However if one of the two packets is stronger by six dB, it will survive."
              
            % https://lora-developers.semtech.com/library/tech-papers-and-guides/lora-and-lorawan/
            % 6 dB  = 3.9811 (approx 4)
            
                  if  pr_h_g_D  < (pr_h_g_I*10^(6/10))
                    discarded = discarded + 1; %% Destructive collision
                  end
                 
               end
            end
            P_SIR(count,n, c) = 1- (discarded/MonteCarlo); %% Success in presence of interfering nodes
        end
    end
end

figure
for count=1:length(D_SNR)
    for i=1:GW_num
        D_tmp=reshape(Distance(count, i, :),1,[]);
        PSIR_tmp=reshape(P_SIR(count,i,:),1,[]);
        if count==1
            PSIRSF11_average(i) = mean(PSIR_tmp(1:length(D_tmp(D_tmp~=0))));
            m1 = plot(D_tmp(D_tmp~=0)/1000,PSIR_tmp(1:length(D_tmp(D_tmp~=0))),'rx');
        end
        hold on 
    end
end
grid on
ylabel('Success probability', 'Interpreter', 'Latex');
xlabel('Distance (km)', 'Interpreter', 'Latex');
legend([m1],{'$P_{SIR}$ (SF11)'},'Interpreter', 'Latex');
z = find(~isnan(PSIRSF11_average));
PSIRSF11_system = mean(PSIRSF11_average(z));

clear D_tmp


for count=1:1:1
    for i=1:1:GW_num
        P_s(:,i,count) =  reshape(P_SNR(count,i, :),1,[]).*reshape(P_SIR(count,i, :),1,[]);
    end
end

figure
for count=1:length(D_SNR)
    for i=1:GW_num
        D_tmp=reshape(Distance(count, i, :),1,[]);
        Ps_tmp=reshape(P_s(:,i,count),1,[]);
        if count==1
            PsF11_average(i) = mean(Ps_tmp(1:length(D_tmp(D_tmp~=0))));
            h1 = plot(D_tmp(D_tmp~=0)/1000,Ps_tmp(1:length(D_tmp(D_tmp~=0))),'rx');
        end
        hold on 
    end
end
grid on
ylabel('Success probability', 'Interpreter', 'Latex');
xlabel('Distance (km)', 'Interpreter', 'Latex');
legend([h1],{'$P_{S}$ (SF11)'},'Interpreter', 'Latex');
z = find(~isnan(PsF11_average));
PsSF11_system = mean(PsF11_average(z));
toc
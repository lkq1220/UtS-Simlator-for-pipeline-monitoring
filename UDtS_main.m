clc;
clear; 
close all;

%% Packet and transmissions parameters
tic
Payload = 23;           % Message payload

[Grid, Center_found,D,Idx] = latlon2UTM('bus.csv');  %Coordinates of Pipeline
[Min_distance, Satellite_subpoint]=min(D);
Pipeline = length(Grid);
Sensors_per_Pipeline = 1;
Total_sensors = Pipeline * Sensors_per_Pipeline;

E = 90:-1:10;               %Elevation Angles
R = 6378e3;                % Radius of earth
H = 550e3;                 %Orbital height 
[Satellite_Link_Farms,Ground_distance,Difference_from_GW_Slant_Range] = Simulations_Distance_Points(R,H,Grid,Satellite_subpoint);
[Distance] = Satellite_Geometry(sort(Ground_distance),H);
select_index = [1 9 26 35 41 50 56 61 65 68 70 72 74:1:81];
% select_index=1:1:81;
Distance= Distance(select_index);
% E_angles = [10 20 30 40 50 60 70 80 90];

R_Time = 900;              % Report period
Period = 60.*60.*1000;     % 1 hour in milliseconds
MonteCarlo = 1e3;    % No. of Iterations

%% Gains and Pt are converted into linear form

Pt = 10^(17/10)/1000;      % Transmit Power of LoRa 17 dBm
Freq_Band = 470e6;         % 470 MHz (frequency band China)

Gr=(10.^((22.6)/10));      %22.6: LoRa Gateway
Gt=(10.^((2.15)/10));      %2.15 dBi: End-device
% eta = 2;
% std=0.1;

%paramters for soil characteristics 
Clay_input=3.7;
VWC_input=0.05;%0.132
Depth=0.6;
% [RealSoilDielectric, ImagSoilDielectric] = clc_die(Clay_input, VWC_input, Freq_Band);

%% Simulator Comparsion U-DtS and DtS
% P_snr=zeros(length(Distance),1);
% P_snr_nocity=zeros(length(Distance),1);
% PL_nocity=zeros(length(Distance),1);
% PL_total=zeros(length(Distance),1);
D_SNR = 10.^([-17.5]./10);
% ToA = [41.216 288.768 991.232];
% ToA = [41.216 144.384 577.536]; %Time on Air for 10 bytes in milliseconds
ToA=[823.296];
Channels=48;
% [RealSoilDielectric, ImagSoilDielectric] = clc_die(Clay_input, VWC_input, Freq_Band);
% [PSNR_UDtS(:,1),PSNR_UDtSFSPL(:,1)] = Probability_SNR (Pt,Gt,Gr,D_SNR,Distance,MonteCarlo,eta,std,RealSoilDielectric,ImagSoilDielectric,Depth);
% 
% PSNR_DtS = gwtosatellite(MonteCarlo,Distance);
% 
% figure
% h(1)=plot(Distance/1000,sort(PSNR_DtS,'descend'),'r-','linewidth',2);
% hold on 
% h(2)=plot(Distance/1000,PSNR_UDtS,'b-','linewidth',2);
% grid on
% set(gca,'fontsize',12);
% ylabel('Dense City Probability', 'Interpreter', 'Latex','fontsize',12);
% xlabel('Distance (km)', 'Interpreter', 'Latex','fontsize',12);
% axis([Distance(1)/1000 Distance(end)/1000 0 1]);
% legend('$P_{SNR}$ (DtS)','$P_{SNR}$ (U-DtS Depth=0.5m)','Interpreter', 'Latex','fontsize',12);

Elevation_Angles = 10:10:90;
Elevation_Angles_steps = 10:1:90;
K_factor = [1.24 3.07 3.24 3.6 3.89 5.63 9.77 17.06 25.11];
k = sort(interp1(Elevation_Angles,K_factor,Elevation_Angles_steps),'descend');


%U-DtS comparsion for different SFs
[RealSoilDielectric, ImagSoilDielectric] = clc_die(Clay_input, VWC_input, Freq_Band);

for count=1:1:length(D_SNR)
    [PSNR_UDtS(count,:)] = Probability_SNR(Pt,Gt,Gr,D_SNR(count),Distance,MonteCarlo,RealSoilDielectric,ImagSoilDielectric,Depth,k,select_index);
end

% figure
% plot(Distance/1000,PSNR_UDtS(1,:),'r-','linewidth',2);
% hold on 
% plot(Distance/1000,PSNR_UDtS(2,:),'b-','linewidth',2);
% hold on 
% plot(Distance/1000,PSNR_UDtS(3,:),'k-','linewidth',2);
% hold on 
% grid on

% ylabel('U-DtS Probability', 'Interpreter', 'Latex','fontsize',14);
% xlabel('Distance from user to satellite (km)','Interpreter','Latex','FontSize', 14);
% axis([Distance(1)/1000 Distance(end)/1000 0 1]);
% legend('$P_{SNR-DtS}$ (SF7)','$P_{SNR-DtS}$ (SF9)','$P_{SNR-DtS}$ (SF11)','Interpreter', 'Latex','fontsize',14);
% set(gca,'fontsize',14);

%Calculate PSIR 

for count=1:1:length(D_SNR)
    
    for c=1:length(Distance)  
        [count c]
        discarded = 0;     % To calculate destructive collision at a given distance(as function of elevation angle)
%       Transmissions = (Period/(R_Time*1000))*Node_number(c); 
        Transmissions = (Period/(R_Time*1000))*Total_sensors; 
        for m=1:MonteCarlo
        simultaneous=0;
        %% calculating the number of interfering signals
        
        Transmission_per_channel=round(Transmissions/Channels);  % Transmission per channel

    
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
				  
                  index = randi(length(Grid),simultaneous,1);
                  
                  Interfering_Nodes_X = Satellite_Link_Farms(select_index(c),index,1);
                  Interfering_Nodes_Y = Satellite_Link_Farms(select_index(c),index,2);
                  
                  Interfering_Nodes =[Interfering_Nodes_X', Interfering_Nodes_Y'];
%                    Location_Nodes_Int = sqrt((Interfering_Nodes(:,1)-Grid(Satellite_subpoint,1)).^2 + (Interfering_Nodes(:,2)-Grid(Satellite_subpoint,2)).^2)';
                 Location_Nodes_Int = sqrt((Interfering_Nodes(:,1)-Satellite_Link_Farms(select_index(c),Satellite_subpoint,1)).^2 + (Interfering_Nodes(:,2)-Satellite_Link_Farms(select_index(c),Satellite_subpoint,2)).^2)';
                                        
   
   				  clear dPropogation
	              clear E_dpro
	              clear kC
	              clear E_AngPro
                  dPropogation=zeros(1,simultaneous); 
                
                  E_dpro=zeros(1,simultaneous);
                       
                 % Distance from interfering nodes to satellite
                 for track=1:1:length(Location_Nodes_Int)
                   dPropogation(1,track) = sqrt(H^2 + Location_Nodes_Int(track).^2); 
                 end

                 for Fin_ang=1:length(Location_Nodes_Int)
                   E_dpro(Fin_ang) = (H*((H+2*R)) - dPropogation(Fin_ang).^2)./(2.*dPropogation(Fin_ang).*R);
                 end

                 E_AngPro = round(asind(E_dpro));
                 kC = interp1(Elevation_Angles,K_factor,E_AngPro);
              %% Log Distance Path loss
%               Lo = eta * 10*log10(4*pi*do*(1/wavelength));
%                 clear L_int
%                 clear L_Lin_int
                clear pr_h_g_I_tmp
                clear pr_h_g_I
%                 L_int=zeros(simultaneous);
%                 L_Lin_int=zeros(1,simultaneous);
                
%               L_int  = Lo + 10*eta*log10( Location_Nodes_Int/do) + std*randn(1,1);
                Lu  = City_U2Aloss(RealSoilDielectric, ImagSoilDielectric,Depth,Freq_Band);
                h1=ones(1,MonteCarlo);
                
                for track=1:1:length(Location_Nodes_Int) 
%                     L_int(track,:) =  us2sloss(Pt,Gt,Gr, MonteCarlo,Distance,pointer,Loss_u,k);
%                     L_Lin_int(track,:)=10.^((L_int(track,:))./10);
                    pr_h_g_I_tmp(track)=mean(us2sloss(Pt,Gt,Gr, MonteCarlo,dPropogation(1,track),91-E_AngPro(track),Lu,kC(track)));
                end
                
                %% total interferance = sum of all the interfering signals 
                 pr_h_g_I = sum(pr_h_g_I_tmp);

                %% Rician fading for desired signals                

                %% Received power of desired signal
%                 L_des  = Lo + 10*eta*log10(Distance(c)/do) + std*randn(1,1);
%                 L_des =  us2sloss(MonteCarlo,Distance(c),82-select_index(c),Lu);
%                 L_des =  U2Aloss(RealSoilDielectric, ImagSoilDielectric,Depth,Distance,eta,Freq_Band)+std*randn(1,1);
%                 L_Lin_des=10.^((L_des)./10);
                
                pr_h_g_D = mean(us2sloss(Pt,Gt,Gr, MonteCarlo,Distance(c),select_index(c),Lu,k(select_index(c))));

                %% Data Collisions and Spreading Factor Orthogonality
                
                 %"Two packets with the same spreading factor arriving at the same time... 
                 ...on the same channel might result in a collision. 
                  ...However if one of the two packets is stronger by six dB, it will survive."
                  
                % https://lora-developers.semtech.com/library/tech-papers-and-guides/lora-and-lorawan/
                % 6 dB  = 3.9811 (approx 4)
                
                      if  pr_h_g_D  < (pr_h_g_I*(10^(6/10)))
                        discarded = discarded + 1; %% Destructive collision
                      end
                n=0;
                end
                
        end
     PSIR_UDtS(count, c) = 1- (discarded/MonteCarlo); %% Success in presence of interfering nodes
    end
end

PSF11_UDtS  = PSNR_UDtS(1,:).*PSIR_UDtS(1,:);
% PSF10_UDtS = PSNR_UDtS(2,:).*PSIR_UDtS(2,:);
% PSF12_UDtS = PSNR_UDtS(3,:).*PSIR_UDtS(3,:);

% figure
% plot(Distance/1000,PSIR_UDtS(1,:),'-r','LineWidth',2);
% hold on
% plot(Distance/1000,PSIR_UDtS(2,:),'-b','LineWidth',2);
% hold on
% plot(Distance/1000,PSIR_UDtS(3,:),'-k','LineWidth',2);
% grid on
% ylabel('U-DtS Probability','Interpreter','Latex','FontSize', 14);
% xlabel('Distance from user to satellite (km)','Interpreter','Latex','FontSize', 14);
% axis([Distance(1)/1000 Distance(end)/1000 0 1]);
% legend('$P_{SIR-DtS}$ (SF7)','$P_{SIR-DtS}$ (SF9)','$P_{SIR-DtS}$ (SF11)','Interpreter', 'Latex','fontsize',14,'Location','best');
% set(gca,'fontsize',14);
% 
% figure
% plot(Distance/1000,PSF7_UDtS,'-r','LineWidth',2);
% hold on
% plot(Distance/1000,PSF10_UDtS,'-b','LineWidth',2);
% hold on
% plot(Distance/1000,PSF12_UDtS,'-k','LineWidth',2);
% grid on
% ylabel('U-DtS Probability','Interpreter','Latex','FontSize', 14);
% xlabel('Distance from user to satellite (km)','Interpreter','Latex','FontSize', 14);
% axis([Distance(1)/1000 Distance(end)/1000 0 1]);
% legend('$P_{S-DtS}$ (SF7)','$P_{S-DtS}$ (SF9)','$P_{S-DtS}$ (SF11)','Interpreter', 'Latex','fontsize',14,'Location','best');
% set(gca,'fontsize',14);

toc
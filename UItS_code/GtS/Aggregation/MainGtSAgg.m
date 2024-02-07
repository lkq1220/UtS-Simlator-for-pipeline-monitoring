clc;
clear; 
close all;

%% Packet and transmissions parameters
tic
Max_MAC_payload = [250 192 31];           % Message payload
Frame_payload = 10;
Frame_header = 7;
Frame_port = 1;
MAC_payload = Frame_payload+Frame_header+Frame_port;
MAC_header = 1;
MAC_MIC=4;
Reduce_packet = floor(Max_MAC_payload/MAC_payload);
PHY_payload = (MAC_payload+MAC_header)*Reduce_packet+MAC_MIC;

GW_num=1000;
[Grid, Idx, GW_center,D_center,Center_found] = latlon2UTM('bus.csv',GW_num); 
[Min_distance, Satellite_subpoint]=min(D_center);

Total_sensors = length(Grid);

E = 90:-1:10;               %Elevation Angles
R = 6378e3;                % Radius of earth
H = 550e3;                 %Orbital height
[Satellite_Link_Farms,Satellite_Link_Center_found] = Simulations_Distance_Points(R,H,GW_center,Center_found);
Distance = height2range(H,1,E);
select_index=1:1:81;
Distance= Distance(select_index);
E_angles = [10 20 30 40 50 60 70 80 90];
K_Factor = [1.24 3.07 3.24 3.6 3.89 5.63 9.77 17.06 25.11];
Elevation_Angles_steps = 10:1:90;
K_value = sort(interp1(E_angles,K_Factor,Elevation_Angles_steps),'descend');

R_Time = 900;              % Report period
Period = 60.*60.*1000;     % 1 hour in milliseconds
MonteCarlo = 1e3;    % No. of Iterations

%% Gains and Pt are converted into linear form

Pt = 10^(17/10)/1000;      % Transmit Power of LoRa 17 dBm
Freq_Band = 470e6;         % 470 MHz (frequency band China)

Gr=(10.^((22.6)/10));      %22.6: LoRa Gateway
% Gr=(10.^((17)/10));      %22.6: HAP Gateway
Gt=(10.^((2.15)/10));      %2.15 dBi: End-device
% eta = 2;
% std=0.1;

%paramters for soil characteristics 
Clay_input=3.7;
VWC_input=0.132;
Depth=0.6;
% [RealSoilDielectric, ImagSoilDielectric] = clc_die(Clay_input, VWC_input, Freq_Band);

%% Simulator Comparsion U-DtS and DtS
D_SNR = 10.^([-6 -12 -17.5]./10); %SF7 9 11

ToA=[394.496 984.064 823.296]; %251Bytes 194Bytes 23Bytes
Channels=48;

for count=1:1:length(D_SNR)
    [PSNR_DtS(count,:)] = Probability_SNR(Pt,Gt,Gr,D_SNR(count),Distance,MonteCarlo,Freq_Band);
end

figure
plot(Distance/1000,PSNR_DtS(1,:),'r-','linewidth',2);
hold on 
plot(Distance/1000,PSNR_DtS(2,:),'b-','linewidth',2);
hold on 
plot(Distance/1000,PSNR_DtS(3,:),'k-','linewidth',2);
hold on 
grid on

ylabel('GW2Satellite Probability', 'Interpreter', 'Latex','fontsize',14);
xlabel('Distance from GW to satellite (km)', 'Interpreter', 'Latex','fontsize',14);
axis([Distance(1)/1000 Distance(end)/1000 0 1]);
legend('$P_{SNR}$ (SF7)','$P_{SNR}$ (SF9)','$P_{SNR}$ (SF11)','Interpreter', 'Latex','fontsize',14);
set(gca,'fontsize',14);

for count=1:1:length(D_SNR)
    
    for c=1:length(Distance) 
        [count, c]
        discarded = 0;     % To calculate destructive collision at a given distance(as function of elevation angle)
        Transmissions = (Period/(R_Time*1000))*Total_sensors/Reduce_packet(count); 
        for m=1:MonteCarlo
        simultaneous=0;
        %% calculating the number of interfering signals
        
        Transmission_per_channel=round(Transmissions/Channels);  % Transmission per channel

    
        clear TX_stamp
        TX_stamp  = rand(1,round(Transmission_per_channel))*3600;   %Generating transmissions
        
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
				  
                  index = randi(GW_num,simultaneous,1);
                  
                  Interfering_Nodes_X = Satellite_Link_Farms(select_index(c),index,1);
                  Interfering_Nodes_Y = Satellite_Link_Farms(select_index(c),index,2);
                  
                  Interfering_Nodes =[Interfering_Nodes_X', Interfering_Nodes_Y'];
%                    Location_Nodes_Int = sqrt((Interfering_Nodes(:,1)-Grid(Satellite_subpoint,1)).^2 + (Interfering_Nodes(:,2)-Grid(Satellite_subpoint,2)).^2)';
                  Location_Nodes_Int = sqrt((Interfering_Nodes(:,1)-Satellite_Link_Center_found(select_index(c),1,1)).^2 + (Interfering_Nodes(:,2)-Satellite_Link_Center_found(select_index(c),1,2)).^2)';           
   
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
                
                
              %% Log Distance Path loss
                clear L_int
                clear L_Lin_int
                clear pr_h_g_I_tmp
                clear pr_h_g_I
                
                for Fin_ang=1:length(Location_Nodes_Int)
                   E_dpro(Fin_ang) = (H*((H+2*R)) - dPropogation(Fin_ang).^2)./(2.*dPropogation(Fin_ang).*R);
                end

                E_AngPro = round(asind(E_dpro));
                k_inter = interp1(E_angles,K_Factor,E_AngPro);
                
                for track=1:1:length(Location_Nodes_Int) 
                    kC = k_inter(track);
                    mu = sqrt(kC./(2.*(kC+1)));    % Mean 
                    sigma = sqrt(1./(2.*(kC+1)));          % Variance 
                    hr=sigma*randn(1,1)+mu;
                    hi=1j*(sigma*randn(1,1)+mu);
                    h_f=(abs(hr+hi)).^2; %FSPL H_F

                    L = -147.55 + 20*log10(Freq_Band)+20*log10(dPropogation(1,track))+0.1*randn(1,1);
                    L_Lin = 10.^(L./10);
                    pr_h_g_I_tmp (track) = Pt.*Gt.*Gr.*h_f.*(1./L_Lin);
                end
                pr_h_g_I=sum(pr_h_g_I_tmp);
          

                %% Received power of desired signal
                muD = sqrt(K_value(c)./(2.*(K_value(c)+1)));    % Mean 
                sigmaD = sqrt(1./(2.*(K_value(c)+1)));          % Variance 
                hrD=sigmaD*randn(1,1)+muD;
                hiD=1j*(sigmaD*randn(1,1)+muD);
                h_fD=(abs(hrD+hiD)).^2; %FSPL H_F

                LD = -147.55 + 20*log10(Freq_Band)+20*log10(Distance(c))+0.1*randn(1,1);
                L_LinD = 10.^(LD./10);
                pr_h_g_D = Pt.*Gt.*Gr.*h_fD.*(1./L_LinD);

                %% Data Collisions and Spreading Factor Orthogonality
                
                 %"Two packets with the same spreading factor arriving at the same time... 
                 ...on the same channel might result in a collision. 
                  ...However if one of the two packets is stronger by six dB, it will survive."
                  
                % https://lora-developers.semtech.com/library/tech-papers-and-guides/lora-and-lorawan/
                % 6 dB  = 3.9811 (approx 4)
                
                  if  pr_h_g_D  < (pr_h_g_I*(10^(6/10)))
                    discarded = discarded + 1; %% Destructive collision
                  end
                     
               end
        end
     PSIR_DtS(count, c) = 1- (discarded/MonteCarlo); %% Success in presence of interfering nodes
    end
end

PSF7_DtS  = PSNR_DtS(1,:).*PSIR_DtS(1,:);
PSF10_DtS = PSNR_DtS(2,:).*PSIR_DtS(2,:);
PSF12_DtS = PSNR_DtS(3,:).*PSIR_DtS(3,:);

figure
plot(Distance/1000,PSIR_DtS(1,:),'-r','LineWidth',2);
hold on
plot(Distance/1000,PSIR_DtS(2,:),'-b','LineWidth',2);
hold on
plot(Distance/1000,PSIR_DtS(3,:),'-k','LineWidth',2);
grid on
ylabel('GW2Satellite Probability','Interpreter','Latex','FontSize', 14);
xlabel('Distance from GW to satellite (km)','Interpreter','Latex','FontSize', 14);
axis([Distance(1)/1000 Distance(end)/1000 0 1]);
legend('$P_{SIR}$ (SF7)','$P_{SIR}$ (SF9)','$P_{SIR}$ (SF11)','Interpreter', 'Latex','fontsize',14);
set(gca,'fontsize',14);

figure
plot(Distance/1000,PSF7_DtS,'-r','LineWidth',2);
hold on
plot(Distance/1000,PSF10_DtS,'-b','LineWidth',2);
hold on
plot(Distance/1000,PSF12_DtS,'-k','LineWidth',2);
grid on
ylabel('GW2Satellite Probability','Interpreter','Latex','FontSize', 14);
xlabel('Distance from GW to satellite (km)','Interpreter','Latex','FontSize', 14);
axis([Distance(1)/1000 Distance(end)/1000 0 1]);
legend('$P_{S}$ (SF7)','$P_{S}$ (SF9)','$P_{S}$ (SF11)','Interpreter', 'Latex','fontsize',14);
set(gca,'fontsize',14);
toc
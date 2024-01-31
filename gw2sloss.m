function [Dense_FSPL] = gw2sloss(MonteCarlo,Distance,pointer)

tic
%% Communication parameters
% Gr=(10.^((22.6)/10));       %22.6 dBi
% Gt=(10.^((2.15)/10));       %2.15 dBi 
% Ptx = 10^(17/10)/1000;      % Transmit Power of LoRa 14 dBm
Frequency = 470e6;           % in Hz
% MonteCarlo = 1e5;

D_SNR=10.^(-20./10);  %Demodulated SNR of LoRa SF12
Noise_Figure=6;
Band_Width=125e3; %LoRa BW
NoisePower=(10.^((-174+Noise_Figure+10*log10(Band_Width))./10))/1000;

%% Satellite Geometry
Elevation_Angles = 10:10:90;
Elevation_Angles_steps = 10:1:90;
% Orbital_height=550e3;
% Distance = height2range(Orbital_height,1,Elevation_Angles_steps);
%% Distance from user to satellite as function of elevation angle
% [Distance]=Satellite_Geometry(Orbital_height,Elevation_Angles_steps);
% E_angles = [10 20 30 40 50 60 70 80 90];
% K_factor = [1.24 3.07 3.24 3.6 3.89 5.63 9.77 17.06 25.11];
% k = interp1(E_angles,K_factor,Elevation_Angles_steps);


%% Convert elevation angle to slant range
% Distance = height2range(Orbital_height,1,Elevation_Angles_steps);
% Distance = [1.9622e+06,1.8816e+06,1.8059e+06,1.735e+06,1.6685e+06,1.6061e+06,1.5476e+06,1.4927e+06,1.4412e+06,1.3927e+06,1.3472e+06,1.3045e+06,1.2642e+06,1.2263e+06,1.1906e+06,1.1569e+06,1.125e+06,1.095e+06,1.0666e+06,1.0397e+06,1.0143e+06,9.9023e+05,9.674e+05,9.4574e+05,9.2518e+05,9.0566e+05,8.871e+05,8.6946e+05,8.5267e+05,8.367e+05,8.2148e+05,8.0699e+05,7.9317e+05,7.8e+05,7.6743e+05,7.5545e+05,7.44e+05,7.3308e+05,7.2266e+05,7.127e+05,7.0319e+05,6.9411e+05,6.8544e+05,6.7715e+05,6.6924e+05,6.6168e+05,6.5447e+05,6.4758e+05,6.4101e+05,6.3474e+05,6.2876e+05,6.2306e+05,6.1764e+05,6.1247e+05,6.0756e+05,6.0289e+05,5.9846e+05,5.9425e+05,5.9027e+05,5.8651e+05,5.8295e+05,5.796e+05,5.7646e+05,5.735e+05,5.7074e+05,5.6816e+05,5.6577e+05,5.6355e+05,5.6151e+05,5.5965e+05,5.5796e+05,5.5643e+05,5.5507e+05,5.5388e+05,5.5284e+05,5.5197e+05,5.5126e+05,5.5071e+05,5.5031e+05,5.5008e+05,5.5e+05];
%% Line-of-sight (LoS) probability
% Distance=load('Distance.mat').Distance; %Use the default data

Rural_LOS_prob = [78.2, 86.9, 91.9, 92.9, 93.5, 94.0, 94.9, 95.2, 99.8];
Urban_LOS_prob = [24.6, 38.6, 49.3, 61.3, 72.6, 80.5, 91.9, 96.8, 99.2];
Dense_LOS_prob = [28.2, 33.1, 39.8, 46.8, 53.7, 61.2, 73.8, 82.0, 98.1];

Rural_LOS_prob_all_angles=interp1(Elevation_Angles,Rural_LOS_prob ,Elevation_Angles_steps)/100; % For all angles, probability 0-1
Urban_LOS_prob_all_angles=interp1(Elevation_Angles,Urban_LOS_prob ,Elevation_Angles_steps)/100; % For all angles, probability 0-1
Dense_LOS_prob_all_angles=interp1(Elevation_Angles,Dense_LOS_prob ,Elevation_Angles_steps)/100; % For all angles, probability 0-1
% figure(1)
% Lindwidth_value=2;
% plot(Elevation_Angles_steps,Rural_LOS_prob_all_angles,'r-',LineWidth=Lindwidth_value);
% hold on
% plot(Elevation_Angles_steps,Urban_LOS_prob_all_angles,'g-',LineWidth=Lindwidth_value);
% plot(Elevation_Angles_steps,Dense_LOS_prob_all_angles,'b-',LineWidth=Lindwidth_value);
% xlabel('Elevation angle','Interpreter','Latex','FontSize', 12)
% ylabel('LOS probability','Interpreter','Latex','FontSize', 12)
% axis([Elevation_Angles_steps(1) Elevation_Angles_steps(end) 0 1])
% legend('Rural','Urban','Dense Urban','Location','southeast','Interpreter','Latex','FontSize', 12);
% grid on

%% LOS Shadowing for all angles
Rural_LOS_shadow_fading = [1.79, 1.14, 1.14, 0.92, 1.42, 1.56, 0.85, 0.72, 0.72];  %Rural LOS Shadowing
Urban_LOS_shadow_fading = 4*ones(1,length(Rural_LOS_shadow_fading));  %Urban LOS Shadowing
Dense_LOS_shadow_fading = [3.5, 3.4, 2.9, 3.0, 3.1, 2.7, 2.5, 2.3, 1.2];  %Dense Urban LOS Shadowing

Rural_LOS_shadow_all_angles=interp1(Elevation_Angles,Rural_LOS_shadow_fading ,Elevation_Angles_steps); % For all angles
Urban_LOS_shadow_all_angles=interp1(Elevation_Angles,Urban_LOS_shadow_fading ,Elevation_Angles_steps); 
Dense_LOS_shadow_all_angles=interp1(Elevation_Angles,Dense_LOS_shadow_fading ,Elevation_Angles_steps); 
%  zero-mean normal distribution with a standard deviation 
%% https://se.mathworks.com/help/stats/normrnd.html

%% NLOS Shadowing for all angles
Rural_NLOS_shadow_fading = [8.93, 9.08, 8.78, 10.25, 10.56, 10.74, 10.17, 11.52, 11.52];  %Rural NLOS Shadowing 80-11.52 90-11.52
Urban_NLOS_shadow_fading = 6*ones(1,length(Rural_NLOS_shadow_fading)); %Urban NLOS Shadowing
Dense_NLOS_shadow_fading = [15.5, 13.9, 12.4, 11.7, 10.6, 10.5, 10.1, 9.2, 9.2]; %Dense Urban NLOS Shadowing

Rural_NLOS_shadow_all_angles=interp1(Elevation_Angles,Rural_NLOS_shadow_fading ,Elevation_Angles_steps); % For all angles
Urban_NLOS_shadow_all_angles=interp1(Elevation_Angles,Urban_NLOS_shadow_fading ,Elevation_Angles_steps); % For all angles
Dense_NLOS_shadow_all_angles=interp1(Elevation_Angles,Dense_NLOS_shadow_fading ,Elevation_Angles_steps); % For all angles

%% Clutter losses for NLOS
Rural_NLOS_clutter_loss = [19.52, 18.17, 18.42, 18.28, 18.63, 17.68, 16.50, 16.30, 16.30];
%The clutter losses for Urban and Dense Urban are the same.
Urban_NLOS_clutter_loss = [34.3, 30.9, 29.0, 27.7, 26.8, 26.2, 25.8, 25.5, 25.5];
Dense_NLOS_clutter_loss = [34.3, 30.9, 29.0, 27.7, 26.8, 26.2, 25.8, 25.5, 25.5]; 

Rural_NLOS_clutter_loss_all_angles=interp1(Elevation_Angles,Rural_NLOS_clutter_loss,Elevation_Angles_steps); % For all angles
Urban_NLOS_clutter_loss_all_angles=interp1(Elevation_Angles,Urban_NLOS_clutter_loss,Elevation_Angles_steps); % For all angles
Dense_NLOS_clutter_loss_all_angles=interp1(Elevation_Angles,Dense_NLOS_clutter_loss,Elevation_Angles_steps); % For all angles

%% Shadowing fading ploting
% figure(2)
% F_LoS(1)=plot(Elevation_Angles_steps,Rural_LOS_shadow_all_angles,'r-',LineWidth=Lindwidth_value);
% hold on
% F_LoS(2)=plot(Elevation_Angles_steps,Urban_LOS_shadow_all_angles,'g-',LineWidth=Lindwidth_value);
% F_LoS(3)=plot(Elevation_Angles_steps,Dense_LOS_shadow_all_angles,'b-',LineWidth=Lindwidth_value);
% 
% F_NLoS(1)=plot(Elevation_Angles_steps,Rural_NLOS_shadow_all_angles,'r--',LineWidth=Lindwidth_value);
% hold on
% F_NLoS(2)=plot(Elevation_Angles_steps,Urban_NLOS_shadow_all_angles,'g--',LineWidth=Lindwidth_value);
% F_NLoS(3)=plot(Elevation_Angles_steps,Dense_NLOS_shadow_all_angles,'b--',LineWidth=Lindwidth_value);
% xlabel('Elevation angle','Interpreter','Latex','FontSize', 12)
% ylabel('Shadowing standard deviation (dB)','Interpreter','Latex','FontSize', 12)
% grid on
% hold off
% L_LoS=legend(F_LoS(1:3),{'Rural','Urban','Dense Urban'},'Location','northeast','NumColumns',1, 'Interpreter', 'Latex','FontSize',12);
% title(L_LoS,'LoS')
% L_LoS_ax=axes('Position',get(gca,'Position'),'Visible','Off');
% hold on
% 
% L_NLoS=legend(L_LoS_ax,F_NLoS(1:3),{'Rural','Urban','Dense Urban'},'Location','northeast','NumColumns',1, 'Interpreter', 'Latex','FontSize',12);
% title(L_NLoS,'NLoS')
% L_NLoS_ax=axes('Position',get(gca,'Position'),'Visible','Off');


%% 3GPP path loss model for rural areas
% FSPL = 32.45 + 20*log10(carrier frequency) + 20*log10(distance) in GHz

%% path loss calculations featuring LOS and NLOS probabilities
i=pointer;
LOS=0;
NLOS=0;
% mu = sqrt(k(i)./(2.*(k(i)+1)));    % Mean 
% sigma = sqrt(1./(2.*(k(i)+1)));          % Variance 
% hr=sigma.*randn(1,MonteCarlo)+mu;
% hi=1j*(sigma.*randn(1,MonteCarlo)+mu);
% h1=(abs(hr+hi)).^2;
for mCarlo=1:1:MonteCarlo
    Prob=rand;
%         if(Prob<=Rural_LOS_prob_all_angles(i))
%             Rural_FSPL(i,mCarlo) = -147.55 + 20*log10(Frequency)+20*log10(Distance(i))+normrnd(0,Rural_LOS_shadow_all_angles(i));
%             LOS = LOS + 1;
%         else
%             Rural_FSPL(i,mCarlo) = -147.55 + 20*log10(Frequency)+20*log10(Distance(i))+normrnd(0,Rural_NLOS_shadow_all_angles(i))+Rural_NLOS_clutter_loss_all_angles(i);
%             NLOS = NLOS + 1;
%         end
% 
%         if(Prob<=Urban_LOS_prob_all_angles(i))
%             Urban_FSPL(i,mCarlo) = -147.55 + 20*log10(Frequency)+20*log10(Distance(i))+normrnd(0,Urban_LOS_shadow_all_angles(i));
%             %LOS = LOS + 1;
%         else
%             Urban_FSPL(i,mCarlo) = -147.55 + 20*log10(Frequency)+20*log10(Distance(i))+normrnd(0,Urban_NLOS_shadow_all_angles(i))+Urban_NLOS_clutter_loss_all_angles(i);
%             %NLOS = NLOS + 1;
%         end

    if(Prob<=Dense_LOS_prob_all_angles(i))
        Dense_FSPL(1,mCarlo) = -147.55 + 20*log10(Frequency)+20*log10(Distance)+normrnd(0,Dense_LOS_shadow_all_angles(i));
        %LOS = LOS + 1;
    else
        Dense_FSPL(1,mCarlo) = -147.55 + 20*log10(Frequency)+20*log10(Distance)+normrnd(0,Dense_NLOS_shadow_all_angles(i))+Dense_NLOS_clutter_loss_all_angles(i);
        %NLOS = NLOS + 1;
    end
end
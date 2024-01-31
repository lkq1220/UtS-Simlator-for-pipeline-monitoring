function [P_SNR] = Probability_SNR (Pt,Gt,Gr,D_SNR,Distance,MonteCarlo,RealSoilDielectric, ImagSoilDielectric,Depth,k,select_index)
rng(1);
%% LOS : rician fading channel
%% PSNR
Noise_Figure=6;           
Band_Width=125e3; %Band Width 
Freq_Band = 470e6;         % 868 MHz (frequency band Europe)

% mu = sqrt(k./(2.*(k+1)));    % Mean 
% sigma = sqrt(1./(2.*(k+1)));          % Variance 
% 
% hr=sigma*randn(1,MonteCarlo)+mu;
% hi=1j*(sigma*randn(1,MonteCarlo)+mu);

%h1 = (sqrt(0.5)*abs((sigma*randn(length(d1),samples) + mu) + (sigma*1j*randn(length(d1),samples))+mu).^2);

% h1=(abs(hr+hi)).^2;
h1=ones(1,MonteCarlo);

NoisePower=(10.^((-174+Noise_Figure+10*log10(Band_Width))/10))/1000;

std=0.1;      
%% Log Distance Path loss
% LosPro=exp(-0.45*cotd(E_angle));
% gm = gmdistribution((1-LosPro).*12,(1-LosPro).*9+(1-LosPro).*LosPro.*12.^2);
% BuildingLoss = random(gm,MonteCarlo);
% L_return=mean(L);
for pointer=1:length(Distance) 
    
%     L_nocity(pointer,:)  = U2Aloss(RealSoilDielectric, ImagSoilDielectric,Depth,Distance(pointer),eta,Freq_Band)+std*randn(1,MonteCarlo);
%     P_SNR_nocity(pointer)=(sum((Pt.*Gt.*Gr.*h1.*(1./10.^(L_nocity(pointer,:)./10))./NoisePower) >= D_SNR))/MonteCarlo;
    Lu  = City_U2Aloss(RealSoilDielectric, ImagSoilDielectric,Depth,Freq_Band);
    pr = us2sloss(Pt,Gt,Gr,MonteCarlo,Distance(pointer),select_index(pointer),Lu,k(select_index(pointer)));%+std*randn(1,MonteCarlo);
    P_SNR(pointer) = (sum((pr./NoisePower) >= D_SNR))/MonteCarlo;
end





end
% end
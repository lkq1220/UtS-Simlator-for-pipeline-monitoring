function [P_SNR] = Probability_SNR (Pt,Gt,Gr,D_SNR,Distance,MonteCarlo,RealSoilDielectric, ImagSoilDielectric,Depth,k,select_index)
rng(1);
%% LOS : rician fading channel
%% PSNR
Noise_Figure=6;           
Band_Width=125e3; %Band Width 
Freq_Band = 470e6;         % 868 MHz (frequency band Europe)

NoisePower=(10.^((-174+Noise_Figure+10*log10(Band_Width))/10))/1000;
   
for pointer=1:length(Distance) 
    Lu  = City_U2Aloss(RealSoilDielectric, ImagSoilDielectric,Depth,Freq_Band);
    pr = us2sloss(Pt,Gt,Gr,MonteCarlo,Distance(pointer),select_index(pointer),Lu,k(select_index(pointer)));%+std*randn(1,MonteCarlo);
    P_SNR(pointer) = (sum((pr./NoisePower) >= D_SNR))/MonteCarlo;
end

end
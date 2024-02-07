function [P_SNR] = Probability_SNR (Pt,Gt,Gr,D_SNR,Distance,MonteCarlo,Freq_Band)
Noise_Figure=6;           
Band_Width=125e3; %Band Width 

Elevation_Angles = 10:10:90;
Elevation_Angles_steps = 10:1:90;
K_factor = [1.24 3.07 3.24 3.6 3.89 5.63 9.77 17.06 25.11];
k = sort(interp1(Elevation_Angles,K_factor,Elevation_Angles_steps),'descend');

NoisePower=(10.^((-174+Noise_Figure+10*log10(Band_Width))/10))/1000;

std = 0.1;   
%% Log Distance Path loss
for pointer=1:length(Distance) 
    mu = sqrt(k(pointer)./(2.*(k(pointer)+1)));    % Mean 
    sigma = sqrt(1./(2.*(k(pointer)+1)));          % Variance 
    hr=sigma*randn(1,MonteCarlo)+mu;
    hi=1j*(sigma*randn(1,MonteCarlo)+mu);
    h_f=(abs(hr+hi)).^2; %FSPL H_F
   
    L(pointer,:) = -147.55 + 20*log10(Freq_Band)+20*log10(Distance(pointer))+std*randn(1,MonteCarlo);
    L_Lin(pointer,:)=10.^(L(pointer,:)./10);
    pr = Pt.*Gt.*Gr.*h_f.*(1./L_Lin(pointer,:));
    P_SNR(pointer) = (sum((pr./NoisePower) >= D_SNR))/MonteCarlo;
end

end
% end
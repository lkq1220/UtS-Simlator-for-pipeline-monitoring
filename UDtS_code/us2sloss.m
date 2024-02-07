function [Pr] = us2sloss(Pt,Gt,Gr, MonteCarlo,Distance,pointer,Loss_u,k)
%% Communication parameters
Frequency = 470e6;           % in Hz

%% Satellite Geometry
Elevation_Angles = 10:10:90;
Elevation_Angles_steps = 10:1:90;
Dense_LOS_prob = [28.2, 33.1, 39.8, 46.8, 53.7, 61.2, 73.8, 82.0, 98.1];
Dense_LOS_prob_all_angles=sort(interp1(Elevation_Angles,Dense_LOS_prob ,Elevation_Angles_steps)/100,'descend'); % For all angles, probability 0-1

%% path loss calculations featuring LOS and NLOS probabilities
i=pointer;
d_BS = 25;

parfor mCarlo=1:MonteCarlo
    Prob=rand;
    if(Prob<=Dense_LOS_prob_all_angles(i))
        mu = sqrt(k./(2.*(k+1)));    % Mean 
        sigma = sqrt(1./(2.*(k+1)));          % Variance 
        hr=sigma.*randn(1,1)+mu;
        hi=1j*(sigma.*randn(1,1)+mu);
        h1=(abs(hr+hi)).^2;
        Dense_FSPL_LOS= 10.^((Loss_u-147.55 + 20*log10(Frequency)+20*log10(Distance))./10);
        Pr(1,mCarlo) = Pt.*Gt.*Gr.*h1.*(1./Dense_FSPL_LOS);
    else
        mu = sqrt(k./(2.*(k+1)));    % Mean 
        sigma = sqrt(1./(2.*(k+1)));          % Variance 
        hr=sigma.*randn(1,1)+mu;
        hi=1j*(sigma.*randn(1,1)+mu);
        h1=(abs(hr+hi)).^2;

        mu_NLOS = 0;
        sigma_NLOS = sqrt(1/2);          % Variance 
        hr_NLOS=sigma_NLOS.*randn(1,1)+mu_NLOS;
        hi_NLOS=1j*(sigma_NLOS.*randn(1,1)+mu_NLOS);
        hNLOS=(abs(hr_NLOS+hi_NLOS)).^2;

       
        La_nlos = 32.4+30*log10(d_BS)+20*log10(Frequency/1e9)+normrnd(0,7.8);
        Dense_FSPL_NLOS= 10.^((La_nlos+Loss_u-147.55 + 20*log10(Frequency)+20*log10(Distance))./10);

        Pr(1,mCarlo) = Pt.*Gt.*Gr.*h1.*hNLOS.*(1./Dense_FSPL_NLOS);
    end
end
end
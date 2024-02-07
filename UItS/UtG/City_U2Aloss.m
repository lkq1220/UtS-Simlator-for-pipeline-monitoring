function [Pr] = City_U2Aloss(realSoilDielectric, imagSoilDielectric,depth,distance,Freq_Band,MonteCarlo,Pt,Gt,Gr)
    m_realSoilDielectric = realSoilDielectric;
    m_imagSoilDielectric = imagSoilDielectric;
    m_frequecny = Freq_Band;
    EPSILON_0 =  8.854187817 * 10.0^-12;
    MIU_0 = 4 * pi * 10^-7;
    omega = 2 * pi * m_frequecny;
    magEle = (MIU_0 * m_realSoilDielectric * EPSILON_0 / 2.0);
    imagReal = sqrt(1 + (m_imagSoilDielectric / m_realSoilDielectric)^2);
    m_Alpha = omega * sqrt(magEle * (imagReal - 1));
    m_Beta =  omega * sqrt(magEle * (imagReal + 1));
    
    theta_i =asin(1.0./sqrt(m_realSoilDielectric));
    r1 = abs(depth) ./ cos(theta_i);
    Lu = 6.4 + 20.*log10(r1) + 20.*log10(m_Beta) + 8.69.*m_Alpha.*r1;
    
    sigma_NLOS = sqrt(1/2);          % Variance 
    
    k_LOS = 10.^(1.8./10);                   % Urban LOS K-factor    
    mu_LOS = sqrt(k_LOS./(2.*(k_LOS+1)));    % Mean 
    sigma_LOS = sqrt(1./(2.*(k_LOS+1)));     % Variance 
    
    h_bs = 25;
    h_node = 1.5;

    r2 = sqrt(distance.^2+(h_bs-h_node).^2);                       %d_3D
    d_BP = 4*(h_bs-1)*(h_node-1)*m_frequecny/3e8;                  %d_break point 
    pro_LoS = (18/distance)+exp(-distance/63)*(1-18/distance);     %LOS pro in urban scenarios 

    parfor mCarlo=1:MonteCarlo
        pro_value = rand;
        if pro_value<=pro_LoS
            hr_LOS=sigma_LOS*randn(1,1)+mu_LOS;
            hi_LOS=1j*(sigma_LOS*randn(1,1)+mu_LOS);
            h_LOS=(abs(hr_LOS+hi_LOS)).^2; 
            if distance<=d_BP
                PL_LOS = 22*log10(r2)+28+20*log10(m_frequecny/1e9)+normrnd(0,4)+Lu;
            else
                PL_LOS = 40*log10(r2)+28+20*log10(m_frequecny/1e9)-9*log10(d_BP^2+(h_bs-h_node)^2)+normrnd(0,4)+Lu;
            end
            Pr(1,mCarlo) = Pt.*Gt.*Gr.*h_LOS./(10.^(PL_LOS./10));
        else
            hr_NLOS=sigma_NLOS*randn(1,1);
            hi_NLOS=1j*(sigma_NLOS*randn(1,1));
            h_NLOS = (abs(hr_NLOS+hi_NLOS)).^2;
            if distance<=d_BP
                La_los = 22*log10(r2)+28+20*log10(m_frequecny/1e9)+normrnd(0,4);
            else
                La_los = 40*log10(r2)+28+20*log10(m_frequecny/1e9)-9*log10(d_BP^2+(h_bs-h_node)^2)+normrnd(0,4);
            end
            La_nlos = 32.4+30*log10(r2)+20*log10(m_frequecny/1e9)+normrnd(0,7.8);
            PL_NLOS = max(La_los,La_nlos)+Lu;
            Pr(1,mCarlo) = Pt.*Gt.*Gr.*h_NLOS./(10.^(PL_NLOS./10));
        end
    end
end


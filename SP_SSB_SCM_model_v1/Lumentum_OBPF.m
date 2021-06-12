% OBPF_prepare
function Tx_Freq_Data = Lumentum_OBPF(Tx_Freq_Data,Freq_Vector,Freq_shift,ParamControl)
    Mat = csvread('LumentumFilter.csv');

    WL = Mat(:,1)*1e-9;
    Pwr = -Mat(:,2);

    c = 3e8;

    f = c./WL;
 
    
    f_center = c/1541.9e-9;
    
    Freq_Vector_new = Freq_Vector+Freq_shift+f_center;
    Pwr_new = interp1(f,Pwr,Freq_Vector_new); 
    attenuation_dB = Pwr_new-max(Pwr_new);
    
    if ParamControl.LSB_or_RSB == 0
        attenuation_dB = fliplr(attenuation_dB);
    end
    attenuation_dB(isnan(attenuation_dB)) =min(attenuation_dB);
%     figure;
%     plot(Freq_Vector/1e9,attenuation_dB); grid on;
%     axis([-125,125,-40,0]);
    attenuation = 10.^(0.05*attenuation_dB); % for field, not for power 20*log10()

    
    Tx_Freq_Data = Tx_Freq_Data.*attenuation;
       
end



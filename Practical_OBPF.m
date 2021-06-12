% OBPF_prepare
function Tx_Freq_Data = Practical_OBPF(Tx_Freq_Data,Freq_Vector,Freq_shift)
    filename = 'filter.txt';
    fileID = fopen(filename,'r');
    formatSpec = '%f %f';
    A = textscan(fileID,formatSpec);
    fclose(fileID);
    c = 3e8;
    WL = A{1}*1e-9;

    f = c./WL;
    Pwr = A{2};
    
%     figure;
%     plot(WL,Pwr);
%     
%     figure;
%     plot(f/1e9,Pwr);
    
    f_center = c/1550.9e-9;
    
    Freq_Vector_new = Freq_Vector+Freq_shift+f_center;
    Pwr_new = interp1(f,Pwr,Freq_Vector_new); 
    attenuation_dB = Pwr_new-max(Pwr_new);
    figure;
    plot(Freq_Vector/1e9,attenuation_dB); grid on;
    axis([-44,44,-40,0]);
    attenuation = 10.^(0.05*attenuation_dB); % for field, not for power 20*log10()

    
    Tx_Freq_Data = Tx_Freq_Data.*attenuation;

       
end



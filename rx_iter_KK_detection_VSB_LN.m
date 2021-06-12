% ***************************************************
% achieved using hilbert transform
% ***************************************************
function Rx_Time_Data = rx_iter_KK_detection_VSB_LN(Rx_I,iter,...
    DTime,ParamVSB, ParamFib, ParamSig,ParamControl)
    
    %% prepare parameters
    totalBaudRate = sum(ParamSig.SC_Baud_Rate);
    roll_off = ParamSig.roll_off;
    GuardBand = ParamSig.GuardBand;
    
    Fiber_Length = ParamFib.FiberLength; 
    Beta2 = ParamFib.Beta2;
    
    %% get VFreq
    Number_of_Samples = length(Rx_I);
    MaxFreq = 0.5/DTime;
    DFreq = 2*MaxFreq/(Number_of_Samples-1);

    VFreq = (-1:Number_of_Samples-2)*DFreq;
    VFreq = VFreq-0.5*max(VFreq);
    VOmeg = 2*pi*VFreq;	
    Disper = -(1j/2)*Beta2*VOmeg.^2*Fiber_Length; 

    Freq_min = GuardBand;
    Freq_max = GuardBand+(1+roll_off)*totalBaudRate;


%%  

    switch ParamControl.OBPF_option
        case 1
            Attenuation_dB = zeros(size(VFreq)); 
            Slope_Stop = ParamVSB.Opt_Flt_offset; 
            Slope_Start =Slope_Stop-ParamVSB.Opt_Flt_Suppression_dB/ParamVSB.Opt_Flt_Slope_dBper10GHz*10e9;
            if min(VFreq) < Slope_Start
                Attenuation_dB(VFreq<=Slope_Start)= -ParamVSB.Opt_Flt_Suppression_dB;
                length_slope = length(Attenuation_dB(VFreq>Slope_Start & VFreq <=Slope_Stop));
                Attenuation_dB(VFreq>Slope_Start & VFreq <=Slope_Stop) = linspace(-ParamVSB.Opt_Flt_Suppression_dB,0,length_slope);
            else
                ParamVSB.Opt_Flt_Suppression_dB = (Slope_Stop-min(VFreq))/10e9*ParamVSB.Opt_Flt_Slope_dBper10GHz;
                length_slope = length(Attenuation_dB(VFreq <=Slope_Stop));
                Attenuation_dB(VFreq <=Slope_Stop) = linspace(-ParamVSB.Opt_Flt_Suppression_dB,0,length_slope);
            end
            
            if ParamControl.LSB_or_RSB
                Attenuation_dB = fliplr(Attenuation_dB);
            end
        case 5
            M = csvread('LumentumFilter.csv');
            WL = M(:,1)*1e-9;
            Pwr = -M(:,2);
            c = 3e8;
            f = c./WL;
            f_center = c/1541.9e-9;

            Freq_Vector_new = VFreq+f_center+ParamVSB.Opt_Flt_offset;
            Pwr_new = interp1(f,Pwr,Freq_Vector_new); 
            Attenuation_dB = Pwr_new-max(Pwr_new);
            if ParamControl.LSB_or_RSB == 0
                Attenuation_dB = fliplr(Attenuation_dB);
            end
            Attenuation_dB(isnan(Attenuation_dB)) =min(Attenuation_dB);
    end
    attenuation = 10.^(0.05*Attenuation_dB);


    
    H2 = fliplr(attenuation);


    H2 = H2.*exp(-1*Disper);
    H2(VFreq<=0) = 1;

    H1 = attenuation;
    H1 = H1.*exp(1*Disper);
    H1(VFreq<=0) = 1;

%     figure;
%     plot(10*log10(abs(H1.*(VFreq>Freq_min &VFreq <=Freq_max)))); hold on;
%     plot(10*log10(abs(H2.*(VFreq>Freq_min &VFreq <=Freq_max)))); hold on;
%     plot(10*log10(abs((H1+H2).*(VFreq>Freq_min &VFreq <=Freq_max))));


    Rx_I0 = Rx_I;

    E_r = sqrt(abs(Rx_I));
    for idx = 1:iter

        S10 = hilbert(real(E_r));


        S1_F = fftshift(fft(S10));   

        
        S1_F = S1_F./(H2+H1);
        S1_F = S1_F.*(VFreq>Freq_min &VFreq <=Freq_max);
  

        S1 = ifft(ifftshift(S1_F.*H1));
        S2 = ifft(ifftshift(S1_F.*H2));
        


%         E_r = sqrt(Rx_I0+1/4*(S1-conj(S1)-S2+conj(S2)).^2);
        E_r = real(sqrt(Rx_I0-(imag(S1)-imag(S2)).^2));
    end
    
    
%     S10 = hilbert(real(E_r));
% 
% 
%     S1_F = fftshift(fft(S10));   
% 
%         
%     S1_F = S1_F./(H2+H1);
%     S1_F = S1_F.*(VFreq>Freq_min &VFreq <=Freq_max);
%     S3 = ifft(ifftshift(S1_F));
%     Rx_Time_Data = S3;
    
    Rx_Time_Data = E_r;
end 
% ***************************************************
% achieved using hilbert transform
% ***************************************************
function [Rx_Time_Data1,Rx_Time_Data2] = rx_iter_KK_detection_twinSSB_LN(Rx_I1,Rx_I2,iter,...
    DTime,ParamVSB, ParamFib, ParamSig,ParamPS)
    
    %% prepare parameters
    totalBaudRate = sum(ParamSig.SC_Baud_Rate);
    roll_off = ParamPS.Alpha;
    GuardBand = ParamSig.GuardBand;
    
    Fiber_Length = ParamFib.FiberLength; 
    Beta2 = ParamFib.Beta2;
    
    %% get VFreq
    Number_of_Samples = length(Rx_I1);
    MaxFreq = 0.5/DTime;
    DFreq = 2*MaxFreq/(Number_of_Samples-1);

    VFreq = (-1:Number_of_Samples-2)*DFreq;
    VFreq = VFreq-0.5*max(VFreq);
    VOmeg = 2*pi*VFreq;	
    Disper = -(1j/2)*Beta2*VOmeg.^2*Fiber_Length; 

    Freq_min = GuardBand;
    Freq_max = GuardBand+(1+roll_off)*totalBaudRate;


%%    
    Freq_Vector = VFreq;
    Attenuation_dB1 = zeros(size(Freq_Vector)); 
    Slope_Stop1 = ParamVSB.Opt_Flt_offset1; 
    Slope_Start1 =Slope_Stop1-ParamVSB.Opt_Flt_Suppression_dB1/ParamVSB.Opt_Flt_Slope_dBper10GHz1*10e9;
    if min(Freq_Vector) < Slope_Start1
        Attenuation_dB1(Freq_Vector<=Slope_Start1)= -ParamVSB.Opt_Flt_Suppression_dB1;
        length_slope1 = length(Attenuation_dB1(Freq_Vector>Slope_Start1 & Freq_Vector <=Slope_Stop1));
        Attenuation_dB1(Freq_Vector>Slope_Start1 & Freq_Vector <=Slope_Stop1) = linspace(-ParamVSB.Opt_Flt_Suppression_dB1,0,length_slope1);
    else
        ParamVSB.Opt_Flt_Suppression_dB1 = -min(Freq_Vector)/10e9*ParamVSB.Opt_Flt_Slope_dBper10GHz1;
        length_slope1 = length(Attenuation_dB1(Freq_Vector <=Slope_Stop1));
        Attenuation_dB1(Freq_Vector <=Slope_Stop1) = linspace(-ParamVSB.Opt_Flt_Suppression_dB1,0,length_slope1);
    end

    attenuation1 = 10.^(0.05*Attenuation_dB1);
    
    
    
    
    Attenuation_dB2 = zeros(size(Freq_Vector)); 
    Slope_Stop2 = ParamVSB.Opt_Flt_offset2; 
    Slope_Start2 =Slope_Stop2-ParamVSB.Opt_Flt_Suppression_dB2/ParamVSB.Opt_Flt_Slope_dBper10GHz2*10e9;
    if min(Freq_Vector) < Slope_Start2
        Attenuation_dB2(Freq_Vector<=Slope_Start2)= -ParamVSB.Opt_Flt_Suppression_dB2;
        length_slope2 = length(Attenuation_dB2(Freq_Vector>Slope_Start2 & Freq_Vector <=Slope_Stop2));
        Attenuation_dB2(Freq_Vector>Slope_Start2 & Freq_Vector <=Slope_Stop2) = linspace(-ParamVSB.Opt_Flt_Suppression_dB2,0,length_slope2);
    else
        ParamVSB.Opt_Flt_Suppression_dB2 = -min(Freq_Vector)/10e9*ParamVSB.Opt_Flt_Slope_dBper10GHz2;
        length_slope2 = length(Attenuation_dB2(Freq_Vector <=Slope_Stop2));
        Attenuation_dB2(Freq_Vector <=Slope_Stop2) = linspace(-ParamVSB.Opt_Flt_Suppression_dB2,0,length_slope2);
    end

    attenuation2 = 10.^(0.05*Attenuation_dB2);
    
    HF1 = attenuation1;
    HF2 = fliplr(attenuation1); 
    
    HF3 = attenuation2;
    HF4 = fliplr(attenuation2);
    
    Hcd = exp(Disper);
    Hcd_conj = exp(-Disper);
    
    H_RSB = ones(size(HF1));
    H_RSB(1:floor(length(HF1)/2)) = 0;
    
    A = Hcd.*HF1;
    B = Hcd_conj.*HF2;
    C = Hcd.*HF4;
    D = Hcd_conj.*HF3;
    
    denom = B.*C-A.*D;
    
    
    b1 = sqrt(Rx_I1);
    b2 = sqrt(Rx_I2);

    for idx = 1:iter

        Fb1_rsb = fftshift(fft(hilbert(real(b1))));
        Fb2_rsb = fftshift(fft(hilbert(real(b2))));
        
        FS1 = H_RSB.*(B.*Fb2_rsb-D.*Fb1_rsb)./denom;
        FS2 = H_RSB.*(C.*Fb1_rsb-A.*Fb2_rsb)./denom;
        FS1 = FS1.*(VFreq>Freq_min &VFreq <=Freq_max);
        FS2 = FS2.*(VFreq>Freq_min &VFreq <=Freq_max);



        s11 = ifft(ifftshift(FS1.*A));
        s12 = ifft(ifftshift(FS2.*B));
        s21 = ifft(ifftshift(FS1.*C));
        s22 = ifft(ifftshift(FS2.*D));
        
        s10 = ifft(ifftshift(FS1));
        s20 = ifft(ifftshift(FS2));

%         E_r = sqrt(Rx_I0+1/4*(S1-conj(S1)-S2+conj(S2)).^2);
        b1 = sqrt(Rx_I1-(imag(s11)-imag(s12)).^2);
        b2 = sqrt(Rx_I2-(imag(s21)-imag(s22)).^2);
    end

    Rx_Time_Data1 = s10;
    Rx_Time_Data2 = s20;
end 
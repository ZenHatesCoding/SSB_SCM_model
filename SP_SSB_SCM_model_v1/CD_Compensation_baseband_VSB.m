function Output_Samples = CD_Compensation_baseband_VSB(Input_Samples, DTime, Fiber_Length, Beta2,Carrier_Offset,...
                          ParamVSB,ParamControl)
 
    Number_of_Samples = length(Input_Samples);
    MaxFreq = 0.5/DTime;
    DFreq = 2*MaxFreq/(Number_of_Samples-1);
    
    VFreq = (-1:Number_of_Samples-2)*DFreq;
    VFreq = VFreq-0.5*max(VFreq);
    VOmeg = 2*pi*VFreq;	
    
    %Dispersion Parameter
    Disper = -(1j/2)*Beta2*(VOmeg+Carrier_Offset*2*pi).^2*Fiber_Length; 
    % with - dispersion 
    % w/o - disp compensation
    % fft exp(-jwt) phase: jw*t-jbeta*z, engineer convention
  
    %dispersion
    Freq_Samples = fftshift(fft(Input_Samples.').');
    
    
    VFreq2 = (-1:4*Number_of_Samples-2)*DFreq;
    VFreq2 = VFreq2-0.5*max(VFreq2);


    switch ParamControl.OBPF_option
        case 1
            Attenuation_dB = zeros(size(VFreq2)); 
            if ParamControl.LSB_or_RSB
                Slope_Stop = ParamVSB.Opt_Flt_offset+Carrier_Offset; 
            else
                Slope_Stop = ParamVSB.Opt_Flt_offset-Carrier_Offset;
            end
            Slope_Start =Slope_Stop-ParamVSB.Opt_Flt_Suppression_dB/ParamVSB.Opt_Flt_Slope_dBper10GHz*10e9;
            if min(VFreq2) < Slope_Start
                Attenuation_dB(VFreq2<=Slope_Start)= -ParamVSB.Opt_Flt_Suppression_dB;
                length_slope = length(Attenuation_dB(VFreq2>Slope_Start & VFreq2 <=Slope_Stop));
                Attenuation_dB(VFreq2>Slope_Start & VFreq2 <=Slope_Stop) = linspace(-ParamVSB.Opt_Flt_Suppression_dB,0,length_slope);
            else
                ParamVSB.Opt_Flt_Suppression_dB = (Slope_Stop-min(VFreq2))/10e9*ParamVSB.Opt_Flt_Slope_dBper10GHz;
                length_slope = length(Attenuation_dB(VFreq2 <=Slope_Stop));
                Attenuation_dB(VFreq2 <=Slope_Stop) = linspace(-ParamVSB.Opt_Flt_Suppression_dB,0,length_slope);
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
            if ParamControl.LSB_or_RSB
                Freq_Vector_new = VFreq2+f_center+ParamVSB.Opt_Flt_offset+Carrier_Offset;
            else
                Freq_Vector_new = VFreq2+f_center+ParamVSB.Opt_Flt_offset-Carrier_Offset;
            end
            Pwr_new = interp1(f,Pwr,Freq_Vector_new); 
            Attenuation_dB = Pwr_new-max(Pwr_new);
           

            Attenuation_dB(isnan(Attenuation_dB)) =min(Attenuation_dB);
            if ParamControl.LSB_or_RSB == 0
                Attenuation_dB = fliplr(Attenuation_dB);
            end
            
    end
    
    
   
    
    attenuation = 10.^(0.05*Attenuation_dB);
    


    flip_attenuation = fliplr(attenuation);
    
    H2 = flip_attenuation(VFreq2<=(max(VFreq) + 2*Carrier_Offset)...
           & VFreq2>=(min(VFreq)-DFreq + 2*Carrier_Offset));

    H2 = H2.*exp(-1*Disper);
    


    
%     H1 = attenuation(VFreq2<=max(VFreq) & VFreq2>=min(VFreq)-DFreq);
    H1 = attenuation(VFreq2<=max(VFreq) & VFreq2>=min(VFreq));
    H1 = H1.*exp(1*Disper);
    
%     figure; plot(VFreq/1e9, [abs(H2);abs(H1)]);
    
    Output_Freq_Samples = Freq_Samples./(H1+H2);
    Output_Samples = ifft(ifftshift(Output_Freq_Samples).').';
    
end


function Output_Samples = CD_Compensation_baseband_VSB_overlap_save(Input_Samples, DTime, Fiber_Length, Beta2,Carrier_Offset,...
                          ParamVSB,ParamRxDSP,ParamControl)
 

    Nfft = ParamRxDSP.Nfft;
    M = ParamRxDSP.CD_tap;
    L = Nfft-M;
    num_block = ceil(length(Input_Samples)/L);
    

    MaxFreq = 0.5/DTime;
    VFreq = linspace(-MaxFreq,MaxFreq,Nfft);
    VOmeg = 2*pi*VFreq;	
    
    %Dispersion Parameter
    Disper = -(1j/2)*Beta2*(VOmeg+Carrier_Offset*2*pi).^2*Fiber_Length; 
    % with - dispersion 
    % w/o - disp compensation
    % fft exp(-jwt) phase: jw*t-jbeta*z, engineer convention
  
    %dispersion

    DFreq = 2*MaxFreq/(Nfft-1);
    
    VFreq2 = (1:4*Nfft)*DFreq;
    
    VFreq2 = VFreq2-0.5*max(VFreq2);
    
    
    switch ParamControl.OBPF_option
        case 1
            Attenuation_dB = zeros(size(VFreq2)); 
            if ParamControl.LSB_or_RSB
               
                Slope_Stop = -ParamVSB.Opt_Flt_offset-Carrier_Offset; 

                Slope_Start =Slope_Stop+ParamVSB.Opt_Flt_Suppression_dB/ParamVSB.Opt_Flt_Slope_dBper10GHz*10e9;
                if max(VFreq2) > Slope_Start
                    Attenuation_dB(VFreq2>=Slope_Start)= -ParamVSB.Opt_Flt_Suppression_dB;
                    length_slope = length(Attenuation_dB(VFreq2<Slope_Start & VFreq2 >=Slope_Stop));
                    Attenuation_dB(VFreq2<Slope_Start & VFreq2 >=Slope_Stop) = linspace(0,-ParamVSB.Opt_Flt_Suppression_dB,length_slope);
                else
                    ParamVSB.Opt_Flt_Suppression_dB = max(VFreq2)/10e9*ParamVSB.Opt_Flt_Slope_dBper10GHz;
                    length_slope = length(Attenuation_dB(VFreq2 >=Slope_Stop));
                    Attenuation_dB(VFreq2 >=Slope_Stop) = linspace(0,-ParamVSB.Opt_Flt_Suppression_dB,length_slope);
                end
            else
            
                Slope_Stop = ParamVSB.Opt_Flt_offset-Carrier_Offset;

                Slope_Start =Slope_Stop-ParamVSB.Opt_Flt_Suppression_dB/ParamVSB.Opt_Flt_Slope_dBper10GHz*10e9;
                if min(VFreq2) < Slope_Start
                    Attenuation_dB(VFreq2<=Slope_Start)= -ParamVSB.Opt_Flt_Suppression_dB;
                    length_slope = length(Attenuation_dB(VFreq2>Slope_Start & VFreq2 <=Slope_Stop));
                    Attenuation_dB(VFreq2>Slope_Start & VFreq2 <=Slope_Stop) = linspace(-ParamVSB.Opt_Flt_Suppression_dB,0,length_slope);
                else
                    ParamVSB.Opt_Flt_Suppression_dB =(Slope_Stop-min(VFreq2))/10e9*ParamVSB.Opt_Flt_Slope_dBper10GHz;
                    length_slope = length(Attenuation_dB(VFreq2 <=Slope_Stop));
                    Attenuation_dB(VFreq2 <=Slope_Stop) = linspace(-ParamVSB.Opt_Flt_Suppression_dB,0,length_slope);
                end
                       
                
            end
        case 5
            Mat = csvread('LumentumFilter.csv');
            WL = Mat(:,1)*1e-9;
            Pwr = -Mat(:,2);
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
           

            
            if ParamControl.LSB_or_RSB == 0
                Attenuation_dB = fliplr(Attenuation_dB);
            end
            Attenuation_dB(isnan(Attenuation_dB)) =min(Attenuation_dB);

    end
    
%     Attenuation_dB = circshift(Attenuation_dB,-5,2);
    attenuation = 10.^(0.05*Attenuation_dB);
    
    
    flip_attenuation = fliplr(attenuation);

    
    H2 = flip_attenuation(VFreq2<=(max(VFreq) + 2*Carrier_Offset)...
           & VFreq2>=(min(VFreq)-DFreq + 2*Carrier_Offset));
%     H2 = circshift(H2,1,2);
    H2 = H2.*exp(-1*Disper);



    H1 = attenuation(VFreq2<=max(VFreq) & VFreq2>=min(VFreq)-DFreq);
    
    H1 = H1.*exp(1*Disper);
    
%     figure; plot(VFreq/1e9, [abs(H2);abs(H1)]);
%     figure; plot(VFreq2/1e9, [abs(flip_attenuation);abs(attenuation)]);
%     xlim([-30 30]);
    
    Output_Samples = [];
    Input_Samples = [Input_Samples, zeros(1,num_block*L-length(Input_Samples))];
    Input_Samples = [zeros(1,M),Input_Samples];
    
    for i = 1:num_block
        Input_Samples_block = Input_Samples((i-1)*L+1:(i-1)*L+Nfft);
        Freq_Samples_block = fftshift(fft(Input_Samples_block));
        Output_Freq_Samples_block = Freq_Samples_block./(H1+H2);    
        Output_Samples_block = ifft(ifftshift(Output_Freq_Samples_block));
        
        Output_Samples = [Output_Samples,Output_Samples_block(M+1:end)];
        
    end  

end


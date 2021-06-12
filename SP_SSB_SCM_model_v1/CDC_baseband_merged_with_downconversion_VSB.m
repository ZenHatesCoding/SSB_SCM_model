function Output_Samples = CDC_baseband_merged_with_downconversion_VSB(Input_Samples, DTime,Carrier_Offset, ParamControl,...
    ParamFib,ParamSig,ParamRxDSP,ParamMod,ParamPD,ParamVSB)
    alpha = ParamRxDSP.FFT_size_ratio;
    num_circshift = -ParamSig.Ncircshift*alpha;
%     num_circshift = 0;
    if ParamRxDSP.KKoverSamp == 64/28
        Nfft = 1024*alpha;
        num_zero_discard = 128*alpha; 

    elseif ParamRxDSP.KKoverSamp == 72/28
        Nfft = 1152*alpha;
        num_zero_discard = 256*alpha; 

    elseif ParamRxDSP.KKoverSamp == 80/28
        Nfft = 1280*alpha;
        num_zero_discard = 384*alpha;
    elseif ParamRxDSP.KKoverSamp == 88/28
        Nfft = 1408*alpha;
        num_zero_discard = 512*alpha;
    elseif ParamRxDSP.KKoverSamp == 96/28
        Nfft = 1536*alpha;
        num_zero_discard = 640*alpha;  
    end

    roll_off = ParamSig.roll_off;
    num_zero_setting = floor((1-(1+roll_off)/ParamRxDSP.KKoverSamp)*Nfft);

    
    M = ParamRxDSP.CD_tap/(Nfft-num_zero_discard)*Nfft;
    
    
    L = Nfft-M;
    num_fft_block = ceil(length(Input_Samples)/L);
    
    Input_Samples = [Input_Samples, zeros(1,num_fft_block*L-length(Input_Samples))];
    
    Input_Samples = [zeros(1,M),Input_Samples];
    
    
    MaxFreq = 0.5/DTime;
    VFreq = linspace(-MaxFreq,MaxFreq,Nfft-num_zero_discard);
    VOmeg = 2*pi*VFreq;	
    
    
    
    Baud_rate = 1/DTime/2;
    Fiber_Length = ParamFib.FiberLength;
    Beta2 = ParamFib.Beta2_ref;
    

    
    %Dispersion Parameter
    Disper = -(1j/2)*Beta2*(VOmeg+Carrier_Offset*2*pi).^2*Fiber_Length; 
    % with - dispersion 
    % w/o - disp compensation
    % fft exp(-jwt) phase: jw*t-jbeta*z, engineer convention
  
    %dispersion

    DFreq = 2*MaxFreq/(Nfft-1-num_zero_discard);
    
    VFreq2 = (1:4*(Nfft-num_zero_discard))*DFreq;
    
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
%             Freq_Vector_new = VFreq2+f_center+ParamVSB.Opt_Flt_offset-Carrier_Offset;
            Pwr_new = interp1(f,Pwr,Freq_Vector_new); 
            Attenuation_dB = Pwr_new-max(Pwr_new);
           

            
            if ParamControl.LSB_or_RSB == 0
                Attenuation_dB = fliplr(Attenuation_dB);
                Attenuation_dB = circshift(Attenuation_dB,-1,2);
            end
            Attenuation_dB(isnan(Attenuation_dB)) =min(Attenuation_dB);

    end
    
%     Attenuation_dB = circshift(Attenuation_dB,-5,2);
    attenuation = 10.^(0.05*Attenuation_dB);
    
    
    flip_attenuation = fliplr(attenuation);

    
    H2 = flip_attenuation(VFreq2<=(max(VFreq) + 2*Carrier_Offset)...
           & VFreq2>=(min(VFreq)-DFreq + 2*Carrier_Offset));


    H2 = H2.*exp(-1*Disper);



    H1 = attenuation(VFreq2<=max(VFreq) & VFreq2>=min(VFreq)-DFreq);

    H1 = H1.*exp(1*Disper);


    %%*********************
    
    Output_Samples = [];
    ParamControl.KK_real_output_or_Not = 0;
    for i = 1:num_fft_block
        Rx_Time_Data_X_block = Input_Samples((i-1)*L+1:(i-1)*L+Nfft);
        if ParamControl.KK_real_output_or_Not == 0
            Rx_Freq_Data_X_block = fftshift(fft(Rx_Time_Data_X_block));
        else
            if i == num_fft_block
                Rx_Time_Data_X_block2 = zeros(1,Nfft);
            else
                Rx_Time_Data_X_block2 = Input_Samples((i)*L+1:(i)*L+Nfft);
            end
            switch ParamControl.KK_real_output_subFFT_num
                case 1
                    Rx_Time_Data_X_block_c = Rx_Time_Data_X_block + 1i*Rx_Time_Data_X_block2;
                case 2
                    Rx_Time_Data_X_block_c = Rx_Time_Data_X_block2 + 1i*Rx_Time_Data_X_block;
            end
            Rx_Freq_Data_X_block_c = fftshift(fft(Rx_Time_Data_X_block_c));
            switch ParamControl.KK_real_output_subFFT_num
                case 1
                    Rx_Freq_Data_X_block_cr = real(Rx_Freq_Data_X_block_c);
                    Rx_Freq_Data_X_block_ci = imag(Rx_Freq_Data_X_block_c);
                    Rx_Freq_Data_X_block_r = 1/2*(Rx_Freq_Data_X_block_cr + [(Rx_Freq_Data_X_block_cr(1)),fliplr((Rx_Freq_Data_X_block_cr(2:end)))]); 
                    Rx_Freq_Data_X_block_i = 1/2*(Rx_Freq_Data_X_block_ci - [(Rx_Freq_Data_X_block_ci(1)),fliplr((Rx_Freq_Data_X_block_ci(2:end)))]); 
                    Rx_Freq_Data_X_block = Rx_Freq_Data_X_block_r + 1i*Rx_Freq_Data_X_block_i;
                case 2
                    Rx_Freq_Data_X_block_cr = real(Rx_Freq_Data_X_block_c);
                    Rx_Freq_Data_X_block_ci = imag(Rx_Freq_Data_X_block_c);
                    Rx_Freq_Data_X_block_r = 1/2*(Rx_Freq_Data_X_block_ci + [(Rx_Freq_Data_X_block_ci(1)),fliplr((Rx_Freq_Data_X_block_ci(2:end)))]); 
                    Rx_Freq_Data_X_block_i = -1/2*(Rx_Freq_Data_X_block_cr - [(Rx_Freq_Data_X_block_cr(1)),fliplr((Rx_Freq_Data_X_block_cr(2:end)))]); 
                    Rx_Freq_Data_X_block = Rx_Freq_Data_X_block_r + 1i*Rx_Freq_Data_X_block_i;
            end
        end

        % down conversion
        Rx_Freq_Data_X_block = circshift(Rx_Freq_Data_X_block,num_circshift,2);
        % LPF
        Rx_Freq_Data_X_block(1:floor(num_zero_setting/2)) = 0;
        Rx_Freq_Data_X_block(end-floor(num_zero_setting/2)+1:end) = 0;
        % resample
        Rx_Freq_Data_X_block(1:num_zero_discard/2) = [];
        Rx_Freq_Data_X_block(end-num_zero_discard/2+1:end) = [];
        
        Freq_Samples_block = Rx_Freq_Data_X_block;
        Output_Freq_Samples_block = Freq_Samples_block./(H1+H2);     
        Output_Samples_block = ifft(ifftshift(Output_Freq_Samples_block));
        
        Output_Samples = [Output_Samples,Output_Samples_block(ParamRxDSP.CD_tap+1:end)];        
    end   

                   

end


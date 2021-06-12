function Output_Samples = CDC_baseband_merged_with_downconversion_VSB_v2(Input_Samples, DTime,Carrier_Offset, ParamControl,...
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

    
    M = ParamRxDSP.CD_tap*Nfft/(Nfft-num_zero_discard);
    
    
    L = Nfft-M;
    num_fft_block = 2*ceil(length(Input_Samples)/L/2);
    
    Input_Samples = [Input_Samples, zeros(1,num_fft_block*L-length(Input_Samples))];
    
    Input_Samples = [zeros(1,M),Input_Samples];
    
    
    MaxFreq = 0.5/DTime;
    K = 1;
    VFreq = linspace(-MaxFreq,MaxFreq,Nfft*K);
    VFreq = VFreq(1:K:end);
    DFreq = VFreq(2)-VFreq(1);
    VOmeg = 2*pi*VFreq;	

  
    
   
    Fiber_Length = ParamFib.FiberLength;
    Beta2 = ParamFib.Beta2_ref;
    

    
    %Dispersion Parameter
    Disper = -(1j/2)*Beta2*(VOmeg+DFreq/16).^2*Fiber_Length;
    
%     Disper(1:Nfft/2) = -Disper(1:Nfft/2);

    % with - dispersion 
    % w/o - disp compensation
    % fft exp(-jwt) phase: jw*t-jbeta*z, engineer convention
  
    %dispersion
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
            Mat = csvread('LumentumFilter.csv');
            WL = Mat(:,1)*1e-9;
            Pwr = -Mat(:,2);
            c = 3e8;
            f = c./WL;
            f_center = c/1541.9e-9;

            Freq_Vector_new = VFreq+f_center+ParamVSB.Opt_Flt_offset;
            Pwr_new = interp1(f,Pwr,Freq_Vector_new); 
            Attenuation_dB = Pwr_new-max(Pwr_new);
            if ParamControl.LSB_or_RSB == 0
                Attenuation_dB = fliplr(Attenuation_dB);
            end
%             Attenuation_dB = fliplr(Attenuation_dB);
            Attenuation_dB(isnan(Attenuation_dB)) =min(Attenuation_dB);
    end

    
    
    %%
%     bw = ParamMod.Modulator_LPF_BW;
%     ord = ParamMod.Modulator_LPF_order;
%     Hmod=exp(-0.5 *log(2)*(VFreq / (bw)).^(2*ord)).';
%     
%     bw = ParamPD.BW_PIN_TIA;
%     ord = ParamPD.BW_PIN_TIA_order;
%     HPD=exp(-0.5 *log(2)*(VFreq / (bw)).^(2*ord)).';
    
    attenuation = 10.^(0.05*Attenuation_dB);

    H2 = fliplr(attenuation);

    

    H2 = H2.*exp(-1*Disper);
%     H2(VFreq<=0) = 1;

    H1 = attenuation;
    H1 = H1.*exp(1*Disper);
%     H1(VFreq<=0) = 1;
    

%     % trick - only for Lumentum filter with fixed filtering response
    if ParamControl.LSB_or_RSB
        load('LSB_h.mat');  
    else
        load('RSB_h.mat');
    end
    H3 = fftshift(fft(bn.',Nfft-num_zero_discard));

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
%         Rx_Freq_Data_X_block = Rx_Freq_Data_X_block./(H1+H2);
%         Rx_Freq_Data_X_block = Rx_Freq_Data_X_block.*exp(1*Disper/3);
        % down conversion
        Rx_Freq_Data_X_block = circshift(Rx_Freq_Data_X_block,num_circshift,2);
        
        
        % LPF
        Rx_Freq_Data_X_block(1:floor(num_zero_setting/2)-0) = 0;
        Rx_Freq_Data_X_block(end-floor(num_zero_setting/2)+1+0:end) = 0;
        % resample
        Rx_Freq_Data_X_block(1:num_zero_discard/2) = [];
        Rx_Freq_Data_X_block(end-num_zero_discard/2+1:end) = [];
        
        % trick
        
        Rx_Freq_Data_X_block = Rx_Freq_Data_X_block.*H3;
        
        Output_Freq_Samples_block = Rx_Freq_Data_X_block;     
        Output_Samples_block = ifft(ifftshift(Output_Freq_Samples_block));
        
        Output_Samples = [Output_Samples,Output_Samples_block(ParamRxDSP.CD_tap+1:end)];        
    end   

%     Output_Samples = conj(Output_Samples);               

end


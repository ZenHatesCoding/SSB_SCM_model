function Output_Samples = CDC_baseband_merged_with_downconversion(Input_Samples, DTime,Carrier_Offset, ParamControl,...
    ParamFib,ParamSig,ParamRxDSP,ParamMod,ParamPD)
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
    Disper = (1j/2)*Beta2*(VOmeg+Carrier_Offset*2*pi).^2*Fiber_Length; 
    % with - dispersion 
    % w/o - disp compensation
    % fft exp(-jwt) phase: jw*t-jbeta*z, engineer convention
    
    VFreq2 = linspace(Carrier_Offset-Baud_rate/2,Carrier_Offset+Baud_rate/2,Nfft/2);
    VFreq3 = linspace(0,100e9,1e4+1);
    bw = ParamMod.Modulator_LPF_BW;
    ord = ParamMod.Modulator_LPF_order;
    H1=exp(-0.5 *log(2)*(VFreq3 / (bw)).^(2*ord)).';
    
    bw = ParamPD.BW_PIN_TIA;
    ord = ParamPD.BW_PIN_TIA_order;
    H2=exp(-0.5 *log(2)*(VFreq3 / (bw)).^(2*ord)).';
    
    H = interp1(VFreq3,H1.*H2,VFreq2);
    H = [ones(1,Nfft/4),H,ones(1,Nfft/4)];
    H = 1;
    
    
    load('SSB_h.mat');
    
    H3 = fftshift(fft(bn.',Nfft-num_zero_discard));
    
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
        
        Output_Freq_Samples_block = Rx_Freq_Data_X_block;
        Output_Freq_Samples_block = Output_Freq_Samples_block.*exp(Disper)./H; 
        
        Output_Freq_Samples_block = Output_Freq_Samples_block.*H3;
        Output_Samples_block = ifft(ifftshift(Output_Freq_Samples_block));
        
        Output_Samples = [Output_Samples,Output_Samples_block(ParamRxDSP.CD_tap+1:end)];        
    end   

                   

end


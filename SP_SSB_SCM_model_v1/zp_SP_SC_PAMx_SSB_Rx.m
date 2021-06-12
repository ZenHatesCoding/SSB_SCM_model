%% By ZP 05/31/2018
% for SCM only
function [BER_X,SNR_X] = zp_SP_SC_PAMx_SSB_Rx(ParamControl,ParamRxDSP,ParamSig,ParamADC,ParamPS,ParamFib,ParamVSB,Rx_I,Tx_Sequence_SC, SC_idx, SE)    
    %% Tool Func
    FT = inline('DTime*fftshift(fft(ifftshift(AtIn)))','AtIn','DTime'); %Analog Fourier transform
    IFT = inline('1/DTime*fftshift(ifft(ifftshift(AfIn)))','AfIn','DTime'); %Inverse Fourier transform
    
    % Carrier_Offset
    SC_Baud_Rate = ParamSig.SC_Baud_Rate;

         
    SC_BW = SC_Baud_Rate*(1+ParamSig.roll_off);
    
    SC_Offset_list = zeros(size(SC_BW));
    SC_Offset_list(1) = SC_BW(1)/2 -sum(SC_BW)/2;
    
    for idx = 2:length(SC_Offset_list)
        SC_Offset_list(idx) = SC_Offset_list(idx-1) + SC_BW(idx-1)/2 + SC_BW(idx)/2;
    end
    
    Carrier_Offset = SC_Offset_list(SC_idx);
    
    % Baud_Rate and rx_Samp_Rate
    Baud_Rate = SC_Baud_Rate(SC_idx);
    rx_Samp_Rate = ParamADC.ADC_Rate;
    
    sum_Baud_Rate = sum(SC_Baud_Rate);
    
    

    %% Time/Freq Vector
    DTime = 1/rx_Samp_Rate;
    Number_of_Samples = length(Rx_I);
    Time_Vector = (0:Number_of_Samples-1)*DTime;
    MaxFreq = 0.5/DTime; 
    DFreq = 2*MaxFreq/(Number_of_Samples-1);
    Freq_Vector = (-(Number_of_Samples-1)/2:(Number_of_Samples-1)/2)'*DFreq;
    
    
    %% Detection
    Rx_I_X = Rx_I(1,:);
    
    

    
    %% KK algorithm
    if ParamControl.Use_DD_or_KK
        Rx_Time_Data_X = Rx_I_X-mean(Rx_I_X);
    else
        
        Rx_I_X = resample(Rx_I_X,sum_Baud_Rate*ParamRxDSP.KKoverSamp,rx_Samp_Rate);
        switch ParamControl.KK_option
            case 0
                Rx_Time_Data_X = rx_KK_detection(Rx_I_X,ParamControl);
            case 1
                Rx_Time_Data_X = rx_KK_woUpSamp_detection(Rx_I_X,ParamControl);
            case 2
                Rx_Time_Data_X = rx_iter_KK_detection(Rx_I_X,ParamRxDSP.numKKiter);
            case 3
                Rx_Time_Data_X = rx_iter_KK_detection_VSB_LN(Rx_I_X,ParamRxDSP.numKKiter,...
                1/(sum_Baud_Rate*ParamRxDSP.KKoverSamp),...
                ParamVSB, ParamFib, ParamSig,ParamPS);
        end

        Rx_Time_Data_X = resample(Rx_Time_Data_X,rx_Samp_Rate,sum_Baud_Rate*ParamRxDSP.KKoverSamp);

        if length(Rx_Time_Data_X)>length(Freq_Vector)
            Rx_Time_Data_X(end+1+length(Freq_Vector)-length(Rx_Time_Data_X):end) = [];
        elseif length(Rx_Time_Data_X)<length(Freq_Vector)
            Rx_Time_Data_X = [Rx_Time_Data_X,zeros(1,length(Freq_Vector)-length(Rx_Time_Data_X))];
        end  
    end
    %% CDC
    if ParamControl.CD_Compensation_or_Not
        Rx_Time_Data_X = CD_Compensation_PAMx_VSB(Rx_Time_Data_X,1/rx_Samp_Rate, ParamVSB,ParamFib,ParamSig);
%         Rx_Time_Data_X = CD_Compensation_baseband_VSB(Rx_Time_Data_X,1/rx_Samp_Rate,ParamFib.FiberLength,ParamFib.Beta2_ref,Carrier_Offset,ParamVSB);    

    end
    Rx_Freq_Data_X = fftshift(fft(Rx_Time_Data_X ));
    Rx_Freq_Data_X(1:floor(length(Rx_Freq_Data_X)/2)) = 0;
    Rx_Time_Data_X = real(ifft(ifftshift(Rx_Freq_Data_X)));
    

    
    
    %% Freq LPF + resample
    if ParamControl.RxDSP_practical_implementation_or_Not
        %% freq down-conversion + LPF + resample
        % only works for single carrier case now
        alpha = ParamRxDSP.FFT_size_ratio;
        switch ParamControl.FEC_option
            case 1
                Nfft = 1024*alpha;                
                num_zero_insertion = 768*alpha; 
                
            case 2
                Nfft = 768*alpha;
                num_zero_insertion = 256*alpha;
            case 3
                Nfft = 768*alpha;
                num_zero_insertion = 256*alpha;

        end
        roll_off = ParamSig.roll_off;
        num_zero_setting = floor((1-(1+roll_off)/(rx_Samp_Rate/Baud_Rate))*Nfft);
        
        num_fft_block = ceil(length(Rx_Time_Data_X)/Nfft);
        Rx_Time_Data_X = [Rx_Time_Data_X,zeros(1,num_fft_block*Nfft-length(Rx_Time_Data_X))];
        Rx_E_Samples_X = [];
        for i = 1:num_fft_block
            Rx_Time_Data_X_block = Rx_Time_Data_X((i-1)*Nfft+1:i*Nfft);
            if ParamControl.FFT_sharing_or_Not == 0
                Rx_Freq_Data_X_block = fftshift(fft(Rx_Time_Data_X_block));
            else
                if i == num_fft_block
                    Rx_Time_Data_X_block2 = zeros(1,Nfft);
                else
                    Rx_Time_Data_X_block2 = Rx_Time_Data_X((i)*Nfft+1:(i+1)*Nfft);
                end
                switch ParamControl.real_output_subFFT_num
                    case 1
                        Rx_Time_Data_X_block_c = Rx_Time_Data_X_block + 1i*Rx_Time_Data_X_block2;
                    case 2
                        Rx_Time_Data_X_block_c = Rx_Time_Data_X_block2 + 1i*Rx_Time_Data_X_block;
                end
                Rx_Freq_Data_X_block_c = fftshift(fft(Rx_Time_Data_X_block_c));
                switch ParamControl.real_output_subFFT_num
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

            % LPF
            Rx_Freq_Data_X_block(1:floor(num_zero_setting/2)) = 0;
            Rx_Freq_Data_X_block(end-floor(num_zero_setting/2)+1:end) = 0;
            % resample
            Rx_Freq_Data_X_block = [zeros(1,num_zero_insertion/2), Rx_Freq_Data_X_block,zeros(1,num_zero_insertion/2)];
            Rx_Time_Data_X_block = ifft(ifftshift(Rx_Freq_Data_X_block));

            Rx_E_Samples_X = [Rx_E_Samples_X, Rx_Time_Data_X_block];
        end   

        Samples_Per_Symbol = 2;
    else
        % freq down conversion
        Rx_Time_Data_X = Rx_Time_Data_X.*exp(-1i*2*pi*Carrier_Offset*Time_Vector);
        % Resample to rx_Samples_Per_Symbol 2Sps
        Samples_Per_Symbol = 2;
        
        ParamPS.OverSamp = Samples_Per_Symbol;
        Rx_Time_Data_X = resample(Rx_Time_Data_X,Baud_Rate*Samples_Per_Symbol,rx_Samp_Rate); 
        DTime = 1/Baud_Rate/Samples_Per_Symbol;
        Number_of_Samples = length(Rx_Time_Data_X);
        Time_Vector = (0:Number_of_Samples-1)*DTime;
        MaxFreq = 0.5/DTime; 
        DFreq = 2*MaxFreq/(Number_of_Samples-1);
        Freq_Vector = (-(Number_of_Samples-1)/2:(Number_of_Samples-1)/2).'*DFreq;


        % Matched filter
        if ParamControl.Pulseshaping_or_Not && ~ParamControl.Use_RC_or_RRC
            Filter_Coef_Rx = Raised_Cosine_Filter(ParamPS);
            Offset = 1+ParamPS.NbSymbols*ParamPS.OverSamp;
            [Rx_E_Samples_X, t2] = rcosflt(Rx_Time_Data_X,ParamPS.Fd,Samples_Per_Symbol,'filter/Fs',Filter_Coef_Rx);
                                  % 'Fs'
                                  % X is input with sample frequency Fs (i.e., the input signal has
                                  % Fs/Fd samples per symbol). In this case the input signal is not
                                  % upsampled from Fd to Fs but is simply filtered by the raised
                                  % cosine filter.
            Rx_E_Samples_X = Rx_E_Samples_X.';
            Rx_E_Samples_X((length(Rx_E_Samples_X)-(Offset-1)+1):end)=[]; 

        else
            Rx_E_Samples_X = Rx_Time_Data_X;
        end
        % LPF to facilitate rrc
        if ParamControl.Pulseshaping_or_Not == 1
            BW_min = -Baud_Rate*(1+ParamSig.roll_off)/2;
            BW_max = Baud_Rate*(1+ParamSig.roll_off)/2;
            BW_seq = round((BW_min-DFreq:DFreq:BW_max+DFreq)./DFreq)+round(Number_of_Samples/2);
            Filter = zeros(size(Rx_E_Samples_X));
            Filter(BW_seq) = 1;
            Rx_F_X = fftshift(fft(Rx_E_Samples_X));
            Rx_E_Samples_X = ifft(ifftshift(Rx_F_X.*Filter));
        end
    end
    
    Rx_Time_Data_X = Rx_E_Samples_X; 


    
    DTime = 1/Baud_Rate/Samples_Per_Symbol;
    Number_of_Samples = length(Rx_Time_Data_X);
    Time_Vector = (0:Number_of_Samples-1)*DTime;
    MaxFreq = 0.5/DTime; 
    DFreq = 2*MaxFreq/(Number_of_Samples-1);
    Freq_Vector = (-(Number_of_Samples-1)/2:(Number_of_Samples-1)/2).'*DFreq;
    
    % up to now single subcarrier at 2 Sps

    
    if ParamControl.Plot_Spectrum_or_Not
        figure;
        
        plot(Freq_Vector/1e9,10*log10(abs(FT(Rx_Time_Data_X,DTime))),'LineWidth',2);
        xlabel('Freq/GHz'); ylabel('Power/dBm');title('Baseband Rx-X @ 2Sps After KK');
    end
    
    Rx_Time_Data_X = pwr_normalization(Rx_Time_Data_X-mean(Rx_Time_Data_X));

  

    
    %% Synchronization
    Tx_Sequence_X = Tx_Sequence_SC(1,:);

    if ParamControl.Synchronization_or_Not
        
        if size(Tx_Sequence_X,2) > 1
            Tx_Sequence_X = Tx_Sequence_X.';
        end
        if size(Rx_E_Samples_X,2) > 1
            Rx_Time_Data_X = Rx_Time_Data_X.';
        end       
        Tx_1sps_X = Tx_Sequence_X;
        Tx_2sps_X = [Tx_1sps_X,Tx_1sps_X].';
        Tx_2sps_X = Tx_2sps_X(:);
%         Tx_2sps_X = upsample(Tx_1sps_X,2,1);

         
        Rx_Time_Data_X = Rx_Time_Data_X(1:length(Tx_2sps_X));
        
        %% cross correlation
        % X-POL
        Tx_Rx_xcorrdata_XX = abs(xcorr(Rx_Time_Data_X,Tx_2sps_X(1:end)));
        if ParamControl.Plot_CrossCorr_or_Not
            figure;plot(Tx_Rx_xcorrdata_XX);        
        end
        [valPeak_XX IndexPeak_XX] = max(Tx_Rx_xcorrdata_XX);   
        SyncError = IndexPeak_XX-length(Tx_2sps_X);
        disp('SyncError:')
        disp(SyncError);

        Rx_Time_Data_X = circshift(Rx_Time_Data_X,-(IndexPeak_XX-length(Tx_2sps_X)));


        % Verifying whether tx and rx are sync
        Tx_Rx_xcorrdata_XX = abs(xcorr(Rx_Time_Data_X,Tx_2sps_X(1:end)));
        if ParamControl.Plot_CrossCorr_or_Not
            figure;plot(Tx_Rx_xcorrdata_XX);
        end
        [valPeak_XX IndexPeak_XX] = max(Tx_Rx_xcorrdata_XX);
        SyncError = IndexPeak_XX-length(Tx_2sps_X);
    else
        Tx_Sequence_X  =  Tx_Sequence_X .'; % trick
    end
    
    Rx_Time_Data_X = pwr_normalization(Rx_Time_Data_X);
    Tx_Sequence_X = pwr_normalization(Tx_Sequence_X);
    
%     Tx_Bits = rx_PAM4_Decode(Tx_Sequence_X);
%     Rx = real(Rx_I_X (1:2:end));
%     Rx = Rx-mean(Rx);
%     Rx_Bits = rx_PAM4_Decode(pwr_normalization(Rx));
%     numError = sum(abs(Rx_Bits-Tx_Bits))  
    
    %% Equalization
    if ParamControl.Equalization_or_Not
        if ParamControl.PLL_or_Not == 0
            M = (ParamRxDSP.LMS_linear_tap-1)/2;
            [bn] = PAMx_LMS_Train(Rx_Time_Data_X(1:ParamRxDSP.Train_Sequence_Length*Samples_Per_Symbol+M),...
                                   Tx_Sequence_X(1:ParamRxDSP.Train_Sequence_Length),...
                                   M, ParamRxDSP.LMS_Train_step, Samples_Per_Symbol,...
                                   ParamControl.Equalization_Plot_Conv_or_Not);
            if ParamControl.Equalization_Plot_Conv_or_Not

                figure;plot(-M:M,abs(bn));title('bn');
                figure;plot(abs(fftshift(fft(bn))));title('bn Freq domain');
                figure;plot(10*log10(abs(fftshift(fft(bn)))));title('bn Freq domain in dB');
            end
            [Rx_Symbols_X] = PAMx_LMS_DD(Rx_Time_Data_X((ParamRxDSP.Train_Sequence_Length*Samples_Per_Symbol+1-M):end),...
                                          bn,SE,... 
                                          ParamRxDSP.LMS_DD_step,ParamControl.Update_Eqtap_or_Not,Samples_Per_Symbol,...
                                          ParamControl.Equalization_Plot_Conv_or_Not);  

            Tx_Symbols_X = Tx_Sequence_X (ParamRxDSP.Train_Sequence_Length+1:end-round(M/Samples_Per_Symbol)); % ignore training symbols   
            Rx_Symbols_X = Rx_Symbols_X(1:length(Tx_Symbols_X));   
        else
            ParamLMS.Trainmu    =  ParamRxDSP.LMS_Train_step;
            ParamLMS.DDmu    = ParamControl.Update_Eqtap_or_Not* ParamRxDSP.LMS_DD_step;
            ParamLMS.NbTapsLMS = ParamRxDSP.LMS_linear_tap;
            ParamLMS.PLLg = 0e-4;
            ParamLMS.PLLg2 = 0e-4;
            ParamLMS.DisplayFigures = ParamControl.Equalization_Plot_Conv_or_Not;
            Rx_Symbols_X = SP_LMS_PLL_PAM4_TSvMOSMANmodified(ParamLMS,Rx_Time_Data_X,Tx_Sequence_X(1:ParamRxDSP.Train_Sequence_Length));
            Rx_Symbols_X = Rx_Symbols_X((ParamRxDSP.Train_Sequence_Length+1):end).';
            M = (ParamRxDSP.LMS_linear_tap-1)/2;
            Tx_Symbols_X = Tx_Sequence_X (ParamRxDSP.Train_Sequence_Length+1:end-round(M/2)); % ignore training symbols   
            Rx_Symbols_X = Rx_Symbols_X(1:length(Tx_Symbols_X));   
        end
        
    else
        Rx_Symbols_X = Rx_Time_Data_X(1+ParamRxDSP.Train_Sequence_Length*Samples_Per_Symbol:Samples_Per_Symbol:end);
        Tx_Symbols_X = Tx_Sequence_X (ParamRxDSP.Train_Sequence_Length+1:end).';
    end
    

    Rx_Symbols_X = (Rx_Symbols_X(ParamRxDSP.Head+1:end-ParamRxDSP.Head)); 
    Rx_Symbols_X = pwr_normalization(Rx_Symbols_X-mean(Rx_Symbols_X));
    % plot rx constellation 
    if ParamControl.Plot_Constellation_or_Not
%         n = ParamSig.Ninterleave;
        n = 1e4;
        if mod(length(Rx_Symbols_X),n)
               Rx_Symbols_X1 = Rx_Symbols_X(1:end-mod(length(Rx_Symbols_X),n));
        else
            Rx_Symbols_X1 = Rx_Symbols_X;
        end
            
        Rx_Symbols_X1 = reshape(Rx_Symbols_X1,n,length(Rx_Symbols_X1)/n);
        Rx_Symbols_X_comp = Rx_Symbols_X1(1:n/2,:)+1i*Rx_Symbols_X1(n/2+1:n,:);
        Rx_Symbols_X_comp = pwr_normalization(Rx_Symbols_X_comp(:).');

        scatterplot(Rx_Symbols_X_comp);
%         scatterplot(pwr_normalization(Rx_Symbols_X(1:1:1e4)+1i*Rx_Symbols_X(1e4+1:1:2e4)));
        axis([-1.5 1.5 -1.5 1.5]);
        title('Constellation');

    end
    Tx_Symbols_X = pwr_normalization(Tx_Symbols_X(ParamRxDSP.Head+1:end-ParamRxDSP.Head)).';  
    
    switch SE
        case 1
            Rx_Bits_X = (Rx_Symbols_X>0);  
            Tx_Bits_X = (Tx_Symbols_X>0);
        case 2
            Rx_Bits_X = rx_PAM4_Decode(Rx_Symbols_X);  
            Tx_Bits_X = rx_PAM4_Decode(Tx_Symbols_X);
        case 2.5
            n = ParamSig.Ninterleave;
            if mod(length(Rx_Symbols_X),n)
                Rx_Symbols_X = Rx_Symbols_X(1:end-mod(length(Rx_Symbols_X),n));
                Tx_Symbols_X = Tx_Symbols_X(1:end-mod(length(Tx_Symbols_X),n));
            end
            
            Rx_Symbols_X = reshape(Rx_Symbols_X,n,length(Rx_Symbols_X)/n);
            Rx_Symbols_X_comp = Rx_Symbols_X(1:n/2,:)+1i*Rx_Symbols_X(n/2+1:n,:);
            Rx_Symbols_X_comp = pwr_normalization(Rx_Symbols_X_comp(:).');
            
            Rx_Bits_X = rx_32QAM_Decode(Rx_Symbols_X_comp); 
            
            Tx_Symbols_X = reshape(Tx_Symbols_X,n,length(Tx_Symbols_X)/n);
            Tx_Symbols_X_comp = Tx_Symbols_X(1:n/2,:)+1i*Tx_Symbols_X(n/2+1:n,:);
            Tx_Symbols_X_comp = pwr_normalization(Tx_Symbols_X_comp(:).');
            
            Tx_Bits_X = rx_32QAM_Decode(Tx_Symbols_X_comp); 
 
        case 3
            Rx_Bits_X = rx_PAM8_Decode(Rx_Symbols_X); 
            Tx_Bits_X = rx_PAM8_Decode(Tx_Symbols_X); 
    end
    
    
    BER_X = sum(abs(Rx_Bits_X-Tx_Bits_X))/length(Rx_Bits_X);
    if ParamControl.Plot_Error_or_Not
        figure; stem(abs(Rx_Bits_X-Tx_Bits_X));
    end
    Var_X = norm(Rx_Symbols_X-Tx_Symbols_X)^2/(length(Rx_Symbols_X));
    SNR_X = 10*log10(1./Var_X);
    if ParamControl.Plot_Noise_or_Not == 1
        noise = Rx_Symbols_X-Tx_Symbols_X;
        figure; histogram(real(noise));
        figure; plot(10*log10(abs(fftshift(fft(noise)))));
    end
    
end
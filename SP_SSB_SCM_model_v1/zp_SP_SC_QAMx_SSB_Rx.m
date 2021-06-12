%% By ZP 05/31/2018
% for SCM only
function [BER_X,SNR_X] = zp_SP_SC_QAMx_SSB_Rx(ParamControl,ParamRxDSP,ParamSig,ParamADC,ParamPS,ParamFib,ParamMod,ParamPD,ParamVSB,Rx_I,Tx_Sequence_SC, SC_idx,SE)    
    %% Tool Func
    FT = inline('DTime*fftshift(fft(ifftshift(AtIn)))','AtIn','DTime'); %Analog Fourier transform
    IFT = inline('1/DTime*fftshift(ifft(ifftshift(AfIn)))','AfIn','DTime'); %Inverse Fourier transform
    
    % Carrier_Offset
    SC_Baud_Rate = ParamSig.SC_Baud_Rate;
    GuardBand = ParamSig.GuardBand;
         
    SC_BW = SC_Baud_Rate*(1+ParamSig.roll_off);
    
    SC_Offset_list = zeros(size(SC_BW));
    SC_Offset_list(1) = SC_BW(1)/2 + GuardBand;
    
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
    
    %% Plot
    if ParamControl.Plot_Spectrum_or_Not
        figure;
        plot(Freq_Vector/1e9,10*log10(abs(FT(Rx_I_X,DTime))),'LineWidth',2);
        xlabel('Freq/GHz'); ylabel('Power/dBm');
        title('Rx-I-X');
    end
    
    if ParamControl.Use_DD_or_KK
        Rx_Time_Data_X = Rx_I_X - mean(Rx_I_X);
    else
        if ParamControl.RxDSP_practical_implementation_or_Not
           if ParamControl.Digital_Resample_Before_KK_or_Not  
                % digital resampling
                if ParamControl.FFT_sharing_or_Not
                    alpha = 2;
                else
                    alpha = 1;
                end

                if ParamControl.FEC_option == 1
                    if ParamRxDSP.KKoverSamp == 72/28
                        Nfft = 1024*alpha;
                        num_zero_insertion = 128*alpha;  
                    elseif ParamRxDSP.KKoverSamp == 80/28
                        Nfft = 1024*alpha;
                        num_zero_insertion = 256*alpha;
                    elseif ParamRxDSP.KKoverSamp == 84/28
                        Nfft = 1024*alpha;
                        num_zero_insertion = 320*alpha; 
                    elseif ParamRxDSP.KKoverSamp == 88/28
                        Nfft = 1024*alpha;
                        num_zero_insertion = 384*alpha;
                    elseif ParamRxDSP.KKoverSamp == 96/28
                        Nfft = 1024*alpha;
                        num_zero_insertion = 512*alpha;
                    elseif ParamRxDSP.KKoverSamp == 104/28
                        Nfft = 1024*alpha;
                        num_zero_insertion = 640*alpha; 
                    elseif ParamRxDSP.KKoverSamp == 112/28
                        Nfft = 1024*alpha;
                        num_zero_insertion = 768*alpha; 
                    end
                else
                    if ParamRxDSP.KKoverSamp == 90/30
                        Nfft = 1280*alpha;
                        num_zero_insertion = 256*alpha;  
                    end
                end
                if (ParamControl.FEC_option == 1 && ParamRxDSP.KKoverSamp ~= 64/28) ||...
                        (ParamControl.FEC_option == 2 && ParamRxDSP.KKoverSamp ~= 75/30)
                    num_fft_block = ceil(length(Rx_I_X)/Nfft);

                    Rx_I_X_len = length(Rx_I_X);
                    Rx_I_X = [Rx_I_X zeros(1,num_fft_block*Nfft-length(Rx_I_X))];
                    Rx_I_X_temp = [];

                    for i = 1:num_fft_block
                        Rx_Time_Data_X_block = Rx_I_X((i-1)*Nfft+1:i*Nfft); 
                        if ParamControl.FFT_sharing_or_Not
                            Rx_Time_Data_X_block = Rx_Time_Data_X_block(1:Nfft/alpha)+1i*Rx_Time_Data_X_block(Nfft/alpha+1:Nfft);
                            Rx_Freq_Data_X_block = fftshift(fft(Rx_Time_Data_X_block));
                            Rx_Freq_Data_X_block = [zeros(1,num_zero_insertion/2/alpha), Rx_Freq_Data_X_block,zeros(1,num_zero_insertion/2/alpha)];
                            Rx_I_X_block = ifft(ifftshift(Rx_Freq_Data_X_block));
                            Rx_I_X_temp = [Rx_I_X_temp, real(Rx_I_X_block),imag(Rx_I_X_block)];
                        else
                            Rx_Freq_Data_X_block = fftshift(fft(Rx_Time_Data_X_block));
                            Rx_Freq_Data_X_block = [zeros(1,num_zero_insertion/2), Rx_Freq_Data_X_block,zeros(1,num_zero_insertion/2)];
                            Rx_I_X_block = ifft(ifftshift(Rx_Freq_Data_X_block));
                            Rx_I_X_temp = [Rx_I_X_temp, Rx_I_X_block];
                        end
                    end
                    Rx_I_X = Rx_I_X_temp; 
                    Rx_I_X = Rx_I_X(1:round(Rx_I_X_len*ParamRxDSP.KKoverSamp/64*28));

                end            
            end
                        
            switch ParamControl.KK_option
                case 0
                    if ParamControl.Hilbert_time_or_Freq ==1
                        Rx_Time_Data_X = rx_KK_detection_FIR_hilbert(Rx_I_X,ParamRxDSP.hilbert_tap+1);
                    else
                        Rx_Time_Data_X = rx_KK_detection_FFT_hilbert(Rx_I_X,ParamRxDSP,ParamControl);
                    end
                case 1
                    Rx_Time_Data_X = rx_upsample_free_KK_detection_FFT_hilbert(Rx_I_X,ParamRxDSP.hilbert_tap,ParamControl);
                case 2
                    Rx_Time_Data_X = rx_iter_KK_detection_FFT_hilbert(Rx_I_X,ParamControl,ParamRxDSP);
                case 3
                    Rx_Time_Data_X = rx_iter_KK_detection_VSB_LN_FFT(Rx_I_X,ParamRxDSP,...
                    1/(sum_Baud_Rate*ParamRxDSP.KKoverSamp),...
                    ParamVSB, ParamFib, ParamSig,ParamControl);
                case 4
                    Rx_Time_Data_X = rx_iter_SSBI_cancel_FFT(Rx_I_X,ParamRxDSP);

            end
%             if ParamControl.KK_option == 0 && ParamControl.CascadeKK_or_Not == 1
%                 Rx_Time_Data_X = rx_iter_KK_detection_v2(Rx_I_X,Rx_Time_Data_X,1);
%             end
            
            
            %% Plot 
            if ParamControl.Plot_Spectrum_or_Not
                figure;
                DTime = 1/ParamRxDSP.KKoverSamp/ParamSig.SC_Baud_Rate;
                MaxFreq = 0.5/DTime; 
                Number_of_Samples = length(Rx_Time_Data_X);
                DFreq = 2*MaxFreq/(Number_of_Samples-1);
                Freq_Vector = (-(Number_of_Samples-1)/2:(Number_of_Samples-1)/2)'*DFreq;

                plot(Freq_Vector/1e9,10*log10(abs(FT(Rx_Time_Data_X,DTime))),'LineWidth',2);
                xlabel('Freq/GHz'); ylabel('Power/dBm');
                title('Rx-Time-Data-X after KK detection');
                ylim([-140 -30]);
            end
        else
            
            if ParamControl.Digital_Resample_Before_KK_or_Not ==1
                Rx_I_X = resample(Rx_I_X,sum_Baud_Rate*ParamRxDSP.KKoverSamp,rx_Samp_Rate);
            end
            
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
                    ParamVSB, ParamFib, ParamSig,ParamControl);
                case 4
                    Rx_Time_Data_X = rx_iter_SSBI_cancel(Rx_I_X,ParamRxDSP.numKKiter);
            end
            
            
            if ParamControl.Digital_Resample_Before_KK_or_Not ==1
                Rx_Time_Data_X = resample(Rx_Time_Data_X,rx_Samp_Rate,sum_Baud_Rate*ParamRxDSP.KKoverSamp);
            end
            if length(Rx_Time_Data_X)>length(Freq_Vector)
                Rx_Time_Data_X(end+1+length(Freq_Vector)-length(Rx_Time_Data_X):end) = [];
            elseif length(Rx_Time_Data_X)<length(Freq_Vector)
                Rx_Time_Data_X = [Rx_Time_Data_X,zeros(1,length(Freq_Vector)-length(Rx_Time_Data_X))];
            end  
            %% Plot 
            if ParamControl.Plot_Spectrum_or_Not
                figure;
                plot(Freq_Vector/1e9,10*log10(abs(FT(Rx_Time_Data_X,DTime))),'LineWidth',2);
                xlabel('Freq/GHz'); ylabel('Power/dBm');
                title('Rx-Time-Data-X after KK detection');
            end
        end
    end
    if ParamControl.RxDSP_practical_implementation_or_Not && ParamControl.MergeCDCwithdownconversion_or_Not

        Samples_Per_Symbol = 2;  
               
        if ParamControl.VSB_or_Not
            DTime = 1/Baud_Rate/ParamRxDSP.KKoverSamp;
            Rx_Time_Data_X = CDC_baseband_merged_with_downconversion_VSB_v2(Rx_Time_Data_X , DTime,Carrier_Offset, ParamControl,...
            ParamFib,ParamSig,ParamRxDSP,ParamMod,ParamPD,ParamVSB);
        else
            DTime = 1/Baud_Rate/2;
            Rx_Time_Data_X = CDC_baseband_merged_with_downconversion(Rx_Time_Data_X , DTime,Carrier_Offset, ParamControl,...
            ParamFib,ParamSig,ParamRxDSP,ParamMod,ParamPD);
        end

    else
        %% Freq down-conversion + LPF + resample
        if ParamControl.RxDSP_practical_implementation_or_Not

            %% freq down-conversion + LPF + resample
            alpha = ParamRxDSP.FFT_size_ratio;
            if ParamControl.VSB_or_Not
    %             num_circshift = -265*alpha;
    %             num_circshift = -252*alpha;
                 num_circshift = -ParamSig.Ncircshift*alpha;
            else
                num_circshift = -ParamSig.Ncircshift*alpha;
            end
            switch ParamControl.FEC_option
                case 1
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
                    elseif ParamRxDSP.KKoverSamp == 104/28
                        Nfft = 1536*alpha;

                        num_zero_discard = 768*alpha;
                    elseif ParamRxDSP.KKoverSamp == 112/28
                        Nfft = 1536*alpha;

                        num_zero_discard = 896*alpha;    
                    end

                    if ParamRxDSP.KKoverSamp == 84/28
                        Nfft = 1536*alpha;
                        num_circshift = -288*alpha;
                        num_zero_discard = 512*alpha;
                    end
                case 2
                    if ParamRxDSP.KKoverSamp == 75/30
                        Nfft = 1280*alpha;
                        num_circshift = -284*alpha;
                        num_zero_discard = 256*alpha;
                    end
                    if ParamRxDSP.KKoverSamp == 90/30
                        Nfft = 1536*alpha;

                        num_circshift = -284*alpha;
                        num_zero_discard = 512*alpha;
                    end
            end
            roll_off = ParamSig.roll_off;
            num_zero_setting = floor((1-(1+roll_off)/ParamRxDSP.KKoverSamp)*Nfft);

            num_fft_block = ceil(length(Rx_Time_Data_X)/Nfft);
            Rx_Time_Data_X = [Rx_Time_Data_X,zeros(1,num_fft_block*Nfft-length(Rx_Time_Data_X))];
            Rx_E_Samples_X = [];
            for i = 1:num_fft_block
                Rx_Time_Data_X_block = Rx_Time_Data_X((i-1)*Nfft+1:i*Nfft);
                if ParamControl.KK_real_output_or_Not == 0
                    Rx_Freq_Data_X_block = fftshift(fft(Rx_Time_Data_X_block));
                else
                    if i == num_fft_block
                        Rx_Time_Data_X_block2 = zeros(1,Nfft);
                    else
                        Rx_Time_Data_X_block2 = Rx_Time_Data_X((i)*Nfft+1:(i+1)*Nfft);
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
            BW_min = -Baud_Rate*(1+ParamSig.roll_off)/2;
            BW_max = Baud_Rate*(1+ParamSig.roll_off)/2;
            BW_seq = round((BW_min-DFreq:DFreq:BW_max+DFreq)./DFreq)+round(Number_of_Samples/2);
            Filter = zeros(size(Rx_E_Samples_X));
            Filter(BW_seq) = 1;
            Rx_F_X = fftshift(fft(Rx_E_Samples_X));
            Rx_E_Samples_X = ifft(ifftshift(Rx_F_X.*Filter));
        end

        Rx_Time_Data_X = Rx_E_Samples_X; 

    %     ParamControl.RxDSP_practical_implementation_or_Not =1;

        %% CD Compensation @ 2Sps
        if ParamControl.CD_Compensation_or_Not
            if ParamControl.RxDSP_practical_implementation_or_Not
                if ParamControl.KK_option ~= 3
                    Rx_Time_Data_X = CD_Compensation_baseband_overlap_save(Rx_Time_Data_X,1/Baud_Rate/2,ParamFib,ParamMod,ParamPD,Carrier_Offset);  
                else
                    Rx_Time_Data_X = CD_Compensation_baseband_VSB_overlap_save(Rx_Time_Data_X,1/Baud_Rate/2,ParamFib.FiberLength,ParamFib.Beta2_ref,Carrier_Offset,ParamVSB,ParamRxDSP,ParamControl);        

                end
            else
                if ParamControl.VSB_or_Not == 1
                    Rx_Time_Data_X = CD_Compensation_baseband_VSB(Rx_Time_Data_X,1/Baud_Rate/2,ParamFib.FiberLength,ParamFib.Beta2_ref,Carrier_Offset,ParamVSB,ParamControl);    
                else
                    Rx_Time_Data_X = CD_Compensation_baseband(Rx_Time_Data_X,1/Baud_Rate/2,ParamFib.FiberLength,ParamFib.Beta2_ref,Carrier_Offset); 
                end

            end
        end
    
    end
    
    DTime = 1/Baud_Rate/Samples_Per_Symbol;
    Number_of_Samples = length(Rx_Time_Data_X);
    Time_Vector = (0:Number_of_Samples-1)*DTime;
    MaxFreq = 0.5/DTime; 
    DFreq = 2*MaxFreq/(Number_of_Samples-1);
    Freq_Vector = (-(Number_of_Samples-1)/2:(Number_of_Samples-1)/2).'*DFreq;

    if ParamControl.Plot_Spectrum_or_Not
        figure;
        
        plot(Freq_Vector/1e9,10*log10(abs(FT(Rx_Time_Data_X,DTime))),'LineWidth',2);
        xlabel('Freq/GHz'); ylabel('Power/dBm');title('Baseband Rx-X @ 2Sps After LPF and CDC');
    end
    
    Rx_Time_Data_X = pwr_normalization(Rx_Time_Data_X-mean(Rx_Time_Data_X));
    
      
    %% Synchronization
    Tx_Sequence_X = Tx_Sequence_SC(1,:);
%     Tx_Sequence_X = conj(Tx_Sequence_SC(1,:));

    
    if ParamControl.Synchronization_or_Not
        if size(Tx_Sequence_X,2) > 1
            Tx_Sequence_X = Tx_Sequence_X.';
        end
        if size(Rx_Time_Data_X,2) > 1
            Rx_Time_Data_X = Rx_Time_Data_X.';
        end       
        Tx_1sps_X = Tx_Sequence_X;
        Tx_2sps_X = upsample(Tx_1sps_X,2,1);
        
        if length(Rx_Time_Data_X)>length(Tx_2sps_X)
            Rx_Time_Data_X = Rx_Time_Data_X(1:length(Tx_2sps_X));
        end
        
        %% cross correlation
        % X-POL
        Tx_Rx_xcorrdata_XX = abs(xcorr(Rx_Time_Data_X,Tx_2sps_X));
        if ParamControl.Plot_CrossCorr_or_Not
            figure;plot(Tx_Rx_xcorrdata_XX);        
        end
        [valPeak_XX IndexPeak_XX] = max(Tx_Rx_xcorrdata_XX);   
        SyncError = IndexPeak_XX-length(Rx_Time_Data_X);
        disp('SyncError:')
        disp(SyncError);

        Rx_Time_Data_X = circshift(Rx_Time_Data_X,-(IndexPeak_XX-length(Rx_Time_Data_X)));

            
        % Verifying whether tx and rx are sync
        Tx_Rx_xcorrdata_XX = abs(xcorr(Rx_Time_Data_X,Tx_2sps_X));
        if ParamControl.Plot_CrossCorr_or_Not
            figure;plot(Tx_Rx_xcorrdata_XX);
        end
        [valPeak_XX IndexPeak_XX] = max(Tx_Rx_xcorrdata_XX);
        SyncError = IndexPeak_XX-length(Rx_Time_Data_X);
        disp('SyncError:')
        disp(SyncError);
        
        Tx_Sequence_X = Tx_Sequence_X.';  
    end
    
    Rx_Time_Data_X = pwr_normalization(Rx_Time_Data_X);
    Tx_Sequence_X = pwr_normalization(Tx_Sequence_X);
    %% Equalization
    if ParamControl.Equalization_or_Not
        if ParamControl.PLL_or_Not == 0
            M = (ParamRxDSP.LMS_linear_tap-1)/2;
            [bn] = QAMx_LMS_Train(Rx_Time_Data_X(1:ParamRxDSP.Train_Sequence_Length*Samples_Per_Symbol+M),...
                                   Tx_Sequence_X(1:ParamRxDSP.Train_Sequence_Length),...
                                   M, ParamRxDSP.LMS_Train_step, Samples_Per_Symbol,...
                                   ParamControl.Equalization_Plot_Conv_or_Not);
            if ParamControl.Equalization_Plot_Conv_or_Not

                figure;plot(-M:M,abs(bn));title('bn');
                figure;plot(abs(fftshift(fft(bn))));title('bn Freq domain');
                figure;plot(10*log10(abs(fftshift(fft(bn)))));title('bn Freq domain in dB');
            end
%             save('SSB_h.mat','bn');
            [Rx_Symbols_X] = QAMx_LMS_DD(Rx_Time_Data_X((ParamRxDSP.Train_Sequence_Length*Samples_Per_Symbol+1-M):end),...
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
            Rx_Symbols_X = SP_LMS_PLL_16QAM_TSvMOSMANmodified(ParamLMS,Rx_Time_Data_X,Tx_Sequence_X(1:ParamRxDSP.Train_Sequence_Length));
            Rx_Symbols_X = Rx_Symbols_X((ParamRxDSP.Train_Sequence_Length+1):end).';
            M = (ParamRxDSP.LMS_linear_tap-1)/2;
            Tx_Symbols_X = Tx_Sequence_X (ParamRxDSP.Train_Sequence_Length+1:end-round(M/2)); % ignore training symbols   
            Rx_Symbols_X = Rx_Symbols_X(1:length(Tx_Symbols_X));   
        end
        
    else
        Rx_Symbols_X = Rx_Time_Data_X(1+ParamRxDSP.Train_Sequence_Length*Samples_Per_Symbol:Samples_Per_Symbol:end);
        Tx_Symbols_X = Tx_Sequence_X (ParamRxDSP.Train_Sequence_Length+1:end);
    end
    

    Rx_Symbols_X = (Rx_Symbols_X(ParamRxDSP.Head+1:end-ParamRxDSP.Head)); 
    Rx_Symbols_X = pwr_normalization(Rx_Symbols_X-mean(Rx_Symbols_X));
    % plot rx constellation 
    if ParamControl.Plot_Constellation_or_Not
        scatterplot(Rx_Symbols_X);
        axis([-1.5 1.5 -1.5 1.5]);
        title('16QAM');
    end
    
    if ParamControl.Add_internal_noise_or_Not
        Rx_Symbols_X = Add_AWGN(Rx_Symbols_X,20,1,1);
    end
    Tx_Symbols_X = pwr_normalization(Tx_Symbols_X(ParamRxDSP.Head+1:end-ParamRxDSP.Head));  
    
    switch SE
        case 5
            Rx_Bits_X = rx_32QAM_Decode(Rx_Symbols_X);  
            Tx_Bits_X = rx_32QAM_Decode(Tx_Symbols_X);
        case 4
            Rx_Bits_X = rx_16QAM_Decode(Rx_Symbols_X);  
            Tx_Bits_X = rx_16QAM_Decode(Tx_Symbols_X);
        case 3
            Rx_Bits_X = rx_8QAM_Decode(Rx_Symbols_X);  
            Tx_Bits_X = rx_8QAM_Decode(Tx_Symbols_X);
        case 2
            Rx_Bits_X = rx_QPSK_Decode(Rx_Symbols_X);  
            Tx_Bits_X = rx_QPSK_Decode(Tx_Symbols_X);
    end
    
    BER_X = sum(abs(Rx_Bits_X-Tx_Bits_X))/length(Rx_Bits_X);
    if ParamControl.Plot_Error_or_Not
        figure; stem(abs(Rx_Bits_X-Tx_Bits_X));
    end
    Var_X = norm(Rx_Symbols_X-Tx_Symbols_X)^2/(length(Rx_Symbols_X));
    SNR_X = 10*log10(1./Var_X);
    
end
function [Tx_Time_Data,Tx_Sequence,ParamPS] = zp_SP_SCM_SSB_Tx(ParamControl,ParamSig,ParamDAC,ParamFib)
    %% Tool Func
    FT = inline('DTime*fftshift(fft(ifftshift(AtIn)))','AtIn','DTime'); %Analog Fourier transform
    IFT = inline('1/DTime*fftshift(ifft(ifftshift(AfIn)))','AfIn','DTime'); %Inverse Fourier transform
    
    %% Overwrite Important Parameters
    DAC_Sample_Rate = ParamDAC.DAC_Rate;
    SC_Baud_Rate = ParamSig.SC_Baud_Rate;
    GuardBand = ParamSig.GuardBand;
    
    num_SC = length(ParamSig.SC);   
    SC_power = ParamSig.SC_power;
    
    SC_BW = SC_Baud_Rate*(1+ParamSig.roll_off);
    SC_Offset = zeros(size(SC_BW));
    SC_Offset(1) = SC_BW(1)/2 + GuardBand;
    for idx = 2:length(SC_Offset)
        SC_Offset(idx) = SC_Offset(idx-1) + SC_BW(idx-1)/2 + SC_BW(idx)/2;
    end
    
    
    %% Dummy Sequence
    [Tx_Time_Data_dummy,Tx_Sequence_dummy,ParamPS_dummy] = zp_SP_QAMx_Tx(ParamControl,ParamDAC, ParamSig,1,2); 
    
    if ParamControl.TxDSP_practical_implementation_or_Not && ParamControl.MergePreEmpwithUpconversion_or_Not
        Tx_Time_Data_dummy = PreEmp_merged_with_upconversion(Tx_Time_Data_dummy,...
                          ParamDAC,ParamSig,ParamControl);
    end
    Tx_Time_Data = zeros(size(Tx_Time_Data_dummy));
    
    DTime = 1/DAC_Sample_Rate;
    Number_of_Samples = length(Tx_Time_Data_dummy);
    Time_Vector = (0:Number_of_Samples-1)*DTime;
    MaxFreq = 0.5/DTime; 
    DFreq = 2*MaxFreq/(Number_of_Samples-1);
    Freq_Vector = (-(Number_of_Samples-1)/2:(Number_of_Samples-1)/2)'*DFreq;
    
    SC_amp_list = sqrt(SC_Baud_Rate/1e9);
    SC_amp_list = SC_amp_list.*SC_power;
    
    

    for idx = 1:num_SC

        [Tx_Time_Data_SC,Tx_Sequence_SC,ParamPS] = zp_SP_QAMx_Tx(ParamControl,ParamDAC, ParamSig,idx ,ParamSig.SC(idx));
        
        if ParamControl.TxDSP_practical_implementation_or_Not
            if ParamControl.MergePreEmpwithUpconversion_or_Not
                Tx_Time_Data_SC = PreEmp_merged_with_upconversion(Tx_Time_Data_SC,...
                          ParamDAC,ParamSig,ParamControl);
%                 Tx_Time_Data_SC = Tx_Time_Data_SC (1:length(Tx_Time_Data_dummy));
            else
                switch ParamControl.FEC_option
                    case 1
                        Nfft = 1024; 
                        num_circshift = ParamSig.Ncircshift;
                    case 2
                        Nfft = 1024; 
                        num_circshift = ParamSig.Ncircshift;
                end
                num_fft_block = ceil(length(Tx_Time_Data_SC)/Nfft);
                Tx_Time_Data_SC = [Tx_Time_Data_SC zeros(1,num_fft_block*Nfft-length(Tx_Time_Data_SC))];
                Tx_Time_Data_SC_temp = [];
                for i = 1:num_fft_block
                    Tx_Time_Data_SC_block = Tx_Time_Data_SC((i-1)*Nfft+1:i*Nfft);
                    Tx_Freq_Data_SC_block = fftshift(fft(Tx_Time_Data_SC_block));
                    Tx_Freq_Data_SC_block = circshift(Tx_Freq_Data_SC_block,num_circshift,2);
                    Tx_Time_Data_SC_block = ifft(ifftshift(Tx_Freq_Data_SC_block));
                    Tx_Time_Data_SC_temp = [Tx_Time_Data_SC_temp, Tx_Time_Data_SC_block];
                end
                Tx_Time_Data_SC = Tx_Time_Data_SC_temp; 
            end
        else
            Tx_Time_Data_SC = Tx_Time_Data_SC.*exp(1i*2*pi*SC_Offset(idx)*Time_Vector);
        end
        
        Tx_Time_Data = Tx_Time_Data + Tx_Time_Data_SC*SC_amp_list(idx);
        Tx_Sequence{idx} = Tx_Sequence_SC;
        
    end
    
    %% pwr_norm
    Tx_Time_Data = pwr_normalization(Tx_Time_Data);
    
    if ParamControl.CD_Pre_Compensation_or_Not
        Tx_Time_Data = CD_Compensation(Tx_Time_Data,DTime,ParamFib.FiberLength,ParamFib.Beta2);
    end
    
    
    if ParamControl.VSB_or_Not
        Tx_Time_Data = real(Tx_Time_Data);
        if ParamControl.Preemphasis_or_Not ...
                && (ParamControl.TxDSP_practical_implementation_or_Not==0||ParamControl.MergePreEmpwithUpconversion_or_Not==0)
            switch ParamControl.DAC_LPF_option
                case 1
                    Tx_Time_Data = real(anti_lpfilt(Tx_Time_Data,ParamDAC.DAC_Rate/1e9,...
                        ParamDAC.DAC_LPF_BW/1e9,ParamDAC.DAC_LPF_type,ParamDAC.DAC_LPF_order));
                   
                case 2
                    Tx_Time_Data = resample(Tx_Time_Data,88e9,ParamDAC.DAC_Rate);
                    Tx_Time_Data = conv(Tx_Time_Data,ParamDAC.Preemphasis_FIR,'same');
                    Tx_Time_Data = resample(Tx_Time_Data,ParamDAC.DAC_Rate,88e9);
            end
            Tx_Time_Data = pwr_normalization(Tx_Time_Data);
            if ParamControl.Plot_Spectrum_or_Not 
                figure;
                plot(Freq_Vector/1e9,max(20*log10(abs(FT(Tx_Time_Data,DTime))),-150),'LineWidth',2);
                xlabel('Freq/GHz'); ylabel('Power/dBm');
                title('DAC waveform after Pre-emphasis');
            end
        end
        Tx_Time_Data = pwr_normalization(Tx_Time_Data);
        if ParamControl.Clipping_or_Not
            Tx_Time_Data = clipFunctionV2(Tx_Time_Data,ParamDAC.clipping_Prob,round(2^ParamDAC.qnbit_DAC));
            Tx_Time_Data = pwr_normalization(Tx_Time_Data);
        end
    else
        if ParamControl.Preemphasis_or_Not...
                && (ParamControl.TxDSP_practical_implementation_or_Not==0||ParamControl.MergePreEmpwithUpconversion_or_Not==0)
            switch ParamControl.DAC_LPF_option
                case 1
                    Tx_Time_Data_I = real(anti_lpfilt(real(Tx_Time_Data),ParamDAC.DAC_Rate/1e9,...
                        ParamDAC.DAC_LPF_BW/1e9,ParamDAC.DAC_LPF_type,ParamDAC.DAC_LPF_order));
                    Tx_Time_Data_Q = real(anti_lpfilt(imag(Tx_Time_Data),ParamDAC.DAC_Rate/1e9,...
                        ParamDAC.DAC_LPF_BW/1e9,ParamDAC.DAC_LPF_type,ParamDAC.DAC_LPF_order));
                case 2
                    Tx_Time_Data = resample(Tx_Time_Data,88e9,ParamDAC.DAC_Rate);
                    Tx_Time_Data_I = conv(real(Tx_Time_Data),ParamDAC.Preemphasis_FIR,'same');
                    Tx_Time_Data_I = resample(Tx_Time_Data_I,ParamDAC.DAC_Rate,88e9);
                    Tx_Time_Data_Q = conv(imag(Tx_Time_Data),ParamDAC.Preemphasis_FIR,'same');
                    Tx_Time_Data_Q = resample(Tx_Time_Data_Q,ParamDAC.DAC_Rate,88e9);

            end
            Tx_Time_Data = Tx_Time_Data_I + 1i*Tx_Time_Data_Q;
            Tx_Time_Data = pwr_normalization(Tx_Time_Data);
            
            if ParamControl.Plot_Spectrum_or_Not
                figure;
                
                plot(Freq_Vector/1e9,max(20*log10(abs(FT(Tx_Time_Data(1:length(Freq_Vector)),DTime))),-150),'LineWidth',2);
                xlabel('Freq/GHz'); ylabel('Power/dBm');
                title('DAC waveform after Preemphasis');
            end
           
        end
        Tx_Time_Data = pwr_normalization(Tx_Time_Data);
        if ParamControl.Clipping_or_Not
            Tx_Time_Data_I = clipFunctionV2(real(Tx_Time_Data),ParamDAC.clipping_Prob,round(2^ParamDAC.qnbit_DAC));
            Tx_Time_Data_Q = clipFunctionV2(imag(Tx_Time_Data),ParamDAC.clipping_Prob,round(2^ParamDAC.qnbit_DAC));
            Tx_Time_Data = pwr_normalization(Tx_Time_Data_I+1i*Tx_Time_Data_Q);
        end
        
    end
    
    %% Plot Tx side Upconverted signal at DAC rate
    if ParamControl.Plot_Spectrum_or_Not
        figure;
        plot(Freq_Vector/1e9,10*log10(abs(FT(Tx_Time_Data,DTime))),'LineWidth',2);
        xlabel('Freq/GHz'); ylabel('Power/dBm');
        title('Signal at DAC Rate');
    end
end
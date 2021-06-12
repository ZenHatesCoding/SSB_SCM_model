%% By ZP 11/02/2019
% implement freq upconvert and resample in Freq domain
% only works for DAC_rate = 64 Gs/s
% Symbol_rate = 28 Gbaud;
%% By ZP 05/30/2018
% Baseband 16 QAM at DAC rate

function [Tx_Time_Data_X,Tx_Sequence_X,ParamPS] = zp_SP_QAMx_Tx(ParamControl,ParamDAC, ParamSig,SC_idx,SE)    
    %% Tool Func
    FT = inline('DTime*fftshift(fft(ifftshift(AtIn)))','AtIn','DTime'); %Analog Fourier transform
    IFT = inline('1/DTime*fftshift(ifft(ifftshift(AfIn)))','AfIn','DTime'); %Inverse Fourier transform

    %% Overwrite Important Parameters
    DAC_Sample_Rate = ParamDAC.DAC_Rate;
    Baud_Rate = ParamSig.SC_Baud_Rate(SC_idx);
    maxSamplesToLoad = ParamDAC.maxSamplesToLoad;
    freq_resolution = ParamDAC.freq_resolution;
    

    %% truncate the symbols according to the baud rate
    if maxSamplesToLoad*Baud_Rate/DAC_Sample_Rate/128 - floor(maxSamplesToLoad*Baud_Rate/DAC_Sample_Rate/128) == 0
        Imax = maxSamplesToLoad * (Baud_Rate/DAC_Sample_Rate);
    else
        Imax = floor(maxSamplesToLoad/128/((DAC_Sample_Rate/freq_resolution) * 1))*128*(DAC_Sample_Rate/freq_resolution * 1) * (Baud_Rate/DAC_Sample_Rate);
    end
    
    Sequence_Length = round(Imax);
    Samples_Per_Symbol = 2; % for pulse shaping
    Number_of_Samples = Sequence_Length*Samples_Per_Symbol;

    %% Time/Freq Vector
    DTime = 1/Baud_Rate/Samples_Per_Symbol;
    Time_Vector = (0:Number_of_Samples-1)*DTime;
    MaxFreq = 0.5/DTime; 
    DFreq = 2*MaxFreq/(Number_of_Samples-1);
    Freq_Vector = (-(Number_of_Samples-1)/2:(Number_of_Samples-1)/2)'*DFreq;

    
    %% Generate Sequence
    M = 2^SE; % 16QAM
    
    Inf_Bits_X = randi([0,1],[1,Sequence_Length*log2(M)]);
    switch SE
        case 5
            Tx_Sequence_X = tx_32QAM_mod(Inf_Bits_X);
        case 4
            Tx_Sequence_X = tx_16QAM_mod(Inf_Bits_X);
        case 3
            Tx_Sequence_X = tx_8QAM_mod(Inf_Bits_X);
        case 2
            Tx_Sequence_X = tx_QPSK_mod(Inf_Bits_X);
    end
    

    %% Pulse shaping                
    if ParamControl.Pulseshaping_or_Not == 1 ...
            && (ParamControl.TxDSP_practical_implementation_or_Not==0||ParamControl.MergePreEmpwithUpconversion_or_Not==0)
        % Pulse shaping Using RRC
        % parameters for pulse shape
        ParamPS.NbSymbols = 32;
        ParamPS.OverSamp = Samples_Per_Symbol; % Fs of rcosflt () 
                                             % should be an integer in this
                                             % simulator
                                             % required sampling rate at Tx
                                             % side for simplicity 
        ParamPS.NbSamples = ParamPS.NbSymbols * ParamPS.OverSamp;
        ParamPS.Alpha = ParamSig.roll_off;  % 0 to 1 rolling factor
        ParamPS.FFTPad= 8;
        ParamPS.Type = ParamControl.Use_RC_or_RRC; 
                       %0 root raised cosine filter; 1 raised cosine filter

        if ParamPS.Alpha == 0
            ParamPS.Alpha = eps;
        end

        ParamPS.Fd = 1;
        Filter_Coef_Tx = Raised_Cosine_Filter(ParamPS);
        Offset = 1+ParamPS.NbSymbols*ParamPS.OverSamp; 

        [PulseShaped_Samples_X, t1] = rcosflt(Tx_Sequence_X,ParamPS.Fd,ParamPS.OverSamp,'filter',Filter_Coef_Tx);

        Tx_Time_Data_X = PulseShaped_Samples_X.';
        Tx_Time_Data_X(1:(Offset-1)) =[]; 
        %  for RRC, both the above line and the following commented 2 lines work
        %  for RC, only the above line works 
        %         Tx_Time_Data_X(1:(Offset-1)/2) =[]; 
        %         Tx_Time_Data_X((length(Tx_Time_Data_X)-(Offset-1)/2+1):end)=[]; 
        


    else 
        % Pulse Shaping NRZ
%         Tx_Time_Data_X = reshape(ones(Samples_Per_Symbol,1)*Tx_Sequence_X,1,Number_of_Samples);
        Tx_Time_Data_X = upsample(Tx_Sequence_X,2,1);
        ParamPS = {};
    end
    

    
    %% resample to DAC_Sample_Rate
    if ParamControl.TxDSP_practical_implementation_or_Not
        if ParamControl.MergePreEmpwithUpconversion_or_Not == 0
            switch ParamControl.FEC_option
                case 1
                    Nfft = 896; 
                    num_zero_insertion = 128;
                case 2
                    Nfft = 1024; 
                    num_zero_insertion = 256;
            end
            num_fft_block = ceil(length(Tx_Time_Data_X)/Nfft);

            Tx_Time_Data_X = [Tx_Time_Data_X zeros(1,num_fft_block*Nfft-length(Tx_Time_Data_X))];
            Tx_Time_Data_X_temp = [];
            for i = 1:num_fft_block
                Tx_Time_Data_X_block = Tx_Time_Data_X((i-1)*Nfft+1:i*Nfft);
                Tx_Freq_Data_X_block = fftshift(fft(Tx_Time_Data_X_block));
                Tx_Freq_Data_X_block = [zeros(1,num_zero_insertion/2), Tx_Freq_Data_X_block,zeros(1,num_zero_insertion/2)];
                Tx_Time_Data_X_block = ifft(ifftshift(Tx_Freq_Data_X_block));
                Tx_Time_Data_X_temp = [Tx_Time_Data_X_temp, Tx_Time_Data_X_block];
            end
            Tx_Time_Data_X = Tx_Time_Data_X_temp; 
        end
    else
        if mod(DAC_Sample_Rate,(1/DTime))
            Tx_Time_Data_X = resample(Tx_Time_Data_X,round(DAC_Sample_Rate),round(1/DTime));
        else
            Tx_Time_Data_X = resample(Tx_Time_Data_X,DAC_Sample_Rate*DTime,1);
        end
    end
    
    DTime = 1/DAC_Sample_Rate;
    Samples_Per_Symbol = DAC_Sample_Rate/Baud_Rate;
    Number_of_Samples = length(Tx_Time_Data_X);
    Time_Vector = (0:Number_of_Samples-1)*DTime;
    MaxFreq = 0.5/DTime; 
    DFreq = 2*MaxFreq/(Number_of_Samples-1);
    Freq_Vector = (-(Number_of_Samples-1)/2:(Number_of_Samples-1)/2)'*DFreq;
    
    
    %% pwr_norm
    Tx_Time_Data_X = pwr_normalization(Tx_Time_Data_X);
        
                                     
end

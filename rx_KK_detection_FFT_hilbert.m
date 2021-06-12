% ***************************************************
% achieved using 32-tap FIR filter
% ***************************************************
function Rx_Time_Data = rx_KK_detection_FFT_hilbert(Rx_I,ParamRxDSP,ParamControl)
    Nfft = ParamRxDSP.Nfft;
    M = ParamRxDSP.hilbert_tap;
    
    sqrt_Rx_I = sqrt(abs(Rx_I));
%     log_sqrt_Rx_I = log(Rx_I);
    log_sqrt_Rx_I = log(sqrt_Rx_I);

    L = Nfft-M;
    num_block = 2*ceil(length(log_sqrt_Rx_I)/L/2);
    
    log_sqrt_Rx_I = [log_sqrt_Rx_I, zeros(1,num_block*L-length(log_sqrt_Rx_I))];

    log_sqrt_Rx_I = [zeros(1,M),log_sqrt_Rx_I];
 
    H = zeros(1,Nfft);
    H(2:Nfft/2) = -1i;
    H(Nfft/2+1:Nfft) = 1i;
    if ParamControl.FFT_sharing_or_Not == 1
        H2 = zeros(1,Nfft);
        H2(2:Nfft/2) = 1;
        H2(Nfft/2+1:Nfft) = -1;
    end
    hilbert_log_sqrt_Rx_I  = [];
    if ParamControl.FFT_sharing_or_Not == 0
        for i = 1:num_block
            Input_Samples_block = log_sqrt_Rx_I((i-1)*L+1:(i-1)*L+Nfft);
            Freq_Samples_block = fft(Input_Samples_block);
            Output_Freq_Samples_block = Freq_Samples_block.*H;
            Output_Samples_block = ifft(Output_Freq_Samples_block);
            hilbert_log_sqrt_Rx_I  = [hilbert_log_sqrt_Rx_I,Output_Samples_block(M+1:end)];

        end
    else
        for i = 1:2:num_block
            Input_Samples_block = log_sqrt_Rx_I((i-1)*L+1:(i-1)*L+Nfft);
            if i == num_block
                Input_Samples_block2 = zeros(1,Nfft);
            else
                Input_Samples_block2 = log_sqrt_Rx_I((i)*L+1:(i)*L+Nfft);
            end
            Freq_Samples_block = fft(Input_Samples_block+1i*Input_Samples_block2);
            Freq_Samples_block_cr = real(Freq_Samples_block);
            Freq_Samples_block_ci = imag(Freq_Samples_block);
            Freq_Samples_block_r = 1/2*(Freq_Samples_block_cr + [(Freq_Samples_block_cr(1)),fliplr((Freq_Samples_block_cr(2:end)))]); 
            Freq_Samples_block_i = 1/2*(Freq_Samples_block_ci - [(Freq_Samples_block_ci(1)),fliplr((Freq_Samples_block_ci(2:end)))]); 
          
            Freq_Samples_block1 =  Freq_Samples_block_r + 1i*Freq_Samples_block_i;
            
            Freq_Samples_block_r = 1/2*(Freq_Samples_block_ci + [(Freq_Samples_block_ci(1)),fliplr((Freq_Samples_block_ci(2:end)))]); 
            Freq_Samples_block_i = -1/2*(Freq_Samples_block_cr - [(Freq_Samples_block_cr(1)),fliplr((Freq_Samples_block_cr(2:end)))]); 

            Freq_Samples_block2 = Freq_Samples_block_r + 1i*Freq_Samples_block_i;
                       
            Output_Freq_Samples_block = (Freq_Samples_block1.*H+Freq_Samples_block2.*H2);
            Output_Samples_block = ifft(Output_Freq_Samples_block);
            hilbert_log_sqrt_Rx_I  = [hilbert_log_sqrt_Rx_I,real(Output_Samples_block(M+1:end)),imag(Output_Samples_block(M+1:end))];
        end
    end
    if length(hilbert_log_sqrt_Rx_I)>length(sqrt_Rx_I)
        hilbert_log_sqrt_Rx_I = hilbert_log_sqrt_Rx_I(1:length(sqrt_Rx_I));

    end
    
    

    if ParamControl.KK_real_output_or_Not == 0
        Rx_Time_Data = sqrt_Rx_I.*exp(1i*hilbert_log_sqrt_Rx_I);
    else
        Rx_Time_Data = sqrt_Rx_I.*cos(hilbert_log_sqrt_Rx_I);
    end
    
    
    
end
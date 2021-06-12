% ***************************************************
% achieved using hilbert transform
% ***************************************************
function Rx_Time_Data = rx_iter_KK_detection_FFT_hilbert(Rx_I,ParamControl,ParamRxDSP)

    Nfft = ParamRxDSP.Nfft;
    M = ParamRxDSP.hilbert_tap;

    L = Nfft-M;
    num_block = ceil(length(Rx_I)/L);
    Rx_I = [Rx_I, zeros(1,num_block*L-length(Rx_I))];

    H = zeros(1,Nfft);
    H(2:Nfft/2) = -1i;
    H(Nfft/2+1:Nfft) = 1i;
    
    if ParamControl.FFT_sharing_or_Not == 1
        H2 = zeros(1,Nfft);
        H2(2:Nfft/2) = 1;
        H2(Nfft/2+1:Nfft) = -1;
    end
    
    E_r = sqrt(Rx_I);
    E_r = [zeros(1,M),E_r];
    for idx = 1:ParamRxDSP.numKKiter        
        E_i = [];
        
        if ParamControl.FFT_sharing_or_Not == 0
            for i = 1:num_block
                Input_Samples_block = E_r((i-1)*L+1:(i-1)*L+Nfft);
                Freq_Samples_block = fft(Input_Samples_block);
                Output_Freq_Samples_block = Freq_Samples_block.*H;
                Output_Samples_block = ifft(Output_Freq_Samples_block);
                E_i  = [E_i,Output_Samples_block(M+1:end)];

            end
        else
            for i = 1:2:num_block
                Input_Samples_block = E_r((i-1)*L+1:(i-1)*L+Nfft);
                if i == num_block
                    Input_Samples_block2 = zeros(1,Nfft);
                else
                    Input_Samples_block2 = E_r((i)*L+1:(i)*L+Nfft);
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

                Output_Freq_Samples_block = 1*(Freq_Samples_block1.*H+Freq_Samples_block2.*H2);
                Output_Samples_block = ifft(Output_Freq_Samples_block);
                E_i  = [E_i,real(Output_Samples_block(M+1:end)),imag(Output_Samples_block(M+1:end))];
            end
        end
        
        if length(E_i)>length(Rx_I)
            E_i = E_i(1:length(Rx_I));
        end
        E_r = sqrt(abs(Rx_I-E_i.^2));
        E_r = [zeros(1,M),E_r];
    end
    
    Rx_Time_Data = E_r(M+1:end);
    
end
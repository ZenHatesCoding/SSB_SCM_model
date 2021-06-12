function y = rx_iter_SSBI_cancel_FFT(x,ParamRxDSP)
    Nfft = ParamRxDSP.Nfft;
    M = ParamRxDSP.hilbert_tap;
    num_iter = ParamRxDSP.numKKiter;
    

    L = Nfft-M+1;
    num_block = ceil(length(x)/L);
    x = [x, zeros(1,num_block*L-length(x))];

    alpha = 3e-3;
    
    
    H = 1*ones(1,Nfft);
    H(1:Nfft/2)=0;

    y = x;

    for idx = 1:num_iter
        x_ssb = [];
        if M > 1
            y = [zeros(1,M-1),y];
        end
        
        for i = 1:num_block
            Input_Samples_block = y((i-1)*L+1:(i-1)*L+Nfft);
            Freq_Samples_block = fftshift(fft(Input_Samples_block));
            Output_Freq_Samples_block = Freq_Samples_block.*H;

            Output_Samples_block = ifft(ifftshift(Output_Freq_Samples_block));
            x_ssb  = [x_ssb,Output_Samples_block(M:end)];

        end
       
        
        
        y = x-alpha*x_ssb.*conj(x_ssb);
        
    end
end
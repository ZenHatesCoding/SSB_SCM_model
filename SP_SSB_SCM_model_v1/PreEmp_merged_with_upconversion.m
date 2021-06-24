function Output_Samples = PreEmp_merged_with_upconversion(Input_Samples,...
                          ParamDAC,ParamSig,ParamControl)
   
    switch ParamControl.FEC_option
        case 1
            Nfft = 896; 
            num_zero_insertion = 128;
            num_circshift = ParamSig.Ncircshift;
        case 2
            Nfft = 960;
            num_zero_insertion = 64;
            num_circshift = ParamSig.Ncircshift;
    end
    
    
    
    
    %% RC filter
    alpha = ParamSig.roll_off;
    Sps = 2;
    Freq_Vector = linspace(-Sps/2,Sps/2,Nfft);
    Hf = zeros(size(Freq_Vector));
    T = 1;
    Hf(Freq_Vector> (alpha-1)/2/T & Freq_Vector < (1-alpha)/2/T) = 1;
    
    Transition_region_1 = Freq_Vector> (1-alpha)/2/T & Freq_Vector < (1+alpha)/2/T;
    Transition_region_2 = Freq_Vector< -(1-alpha)/2/T & Freq_Vector > -(1+alpha)/2/T;
    Hf(Transition_region_1) = (T/2*(1+cos(pi*T/alpha*(Freq_Vector(Transition_region_1)-(1-alpha)/2/T))));
    Hf(Transition_region_2) = (T/2*(1+cos(pi*T/alpha*(-Freq_Vector(Transition_region_2)+(alpha-1)/2/T))));
    
    %% Gaussian filter
    s_rate = ParamDAC.DAC_Rate/1e9;
    Xlen = Nfft+num_zero_insertion;
    f= linspace(0,s_rate/2,Xlen/2);
    bw = ParamDAC.DAC_LPF_BW/1e9;
    ord = ParamDAC.DAC_LPF_order;
    HLP=exp(-0.5 *log(2)*(f / (bw)).^(2*ord)).';
    if mod(Xlen,2) == 0
        HLP = [HLP;flipud(HLP)];
    else
        HLP = [HLP;flipud(HLP(2:end))];
    end
    HLP = HLP.'; % row vector
    HLP = fftshift(HLP);
    if ParamControl.FEC_option == 1
        M = 14;
    else
        M = 15;
    end
    M1 = M*(Nfft+num_zero_insertion)/Nfft;
    L = Nfft-M;
    
    num_fft_block = 2*ceil(length(Input_Samples)/L/2);
    Input_Samples = [Input_Samples, zeros(1,num_fft_block*L-length(Input_Samples))];
    Input_Samples = [zeros(1,M),Input_Samples];
    Output_Samples = [];
    
    
    
    
    for i = 1:num_fft_block
        Tx_Time_Data_X_block = Input_Samples((i-1)*L+1:(i-1)*L+Nfft);
        Tx_Freq_Data_X_block = fftshift(fft(Tx_Time_Data_X_block));
        Tx_Freq_Data_X_block = Tx_Freq_Data_X_block.*Hf;
        
        
        Tx_Freq_Data_X_block = [zeros(1,num_zero_insertion/2), Tx_Freq_Data_X_block,zeros(1,num_zero_insertion/2)];

        Tx_Freq_Data_X_block = circshift(Tx_Freq_Data_X_block,num_circshift,2);
                
        
        Tx_Freq_Data_X_block = Tx_Freq_Data_X_block./HLP;
        Tx_Time_Data_X_block = ifft(ifftshift(Tx_Freq_Data_X_block));
        Output_Samples = [Output_Samples, Tx_Time_Data_X_block(M1+1:end)];
    end
    
    
    
end
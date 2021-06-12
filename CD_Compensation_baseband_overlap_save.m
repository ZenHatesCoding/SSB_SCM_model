function Output_Samples = CD_Compensation_baseband_overlap_save(Input_Samples, DTime, ParamFib,ParamMod,ParamPD,Carrier_Offset)
    Nfft = 896;
    Baud_rate = 1/DTime/2;
    Fiber_Length = ParamFib.FiberLength;
    Beta2 = ParamFib.Beta2_ref;
%     M = ceil(0.032*(Baud_rate/1e9)^2*Fiber_Length/1000*17/1000);
    M = 28;
    L = Nfft-M+1;
    num_block = ceil(length(Input_Samples)/L);
    
    Input_Samples = [Input_Samples, zeros(1,num_block*L-length(Input_Samples))];
    
    Input_Samples = [zeros(1,M-1),Input_Samples];
    
    MaxFreq = 0.5/DTime;
    VFreq = linspace(-MaxFreq,MaxFreq,Nfft);
    VOmeg = 2*pi*VFreq;	
    
    
    
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
    %dispersion
    Output_Samples = [];
    for i = 1:num_block
        Input_Samples_block = Input_Samples((i-1)*L+1:(i-1)*L+Nfft);
        Freq_Samples_block = fftshift(fft(Input_Samples_block));
        Output_Freq_Samples_block = Freq_Samples_block.*exp(Disper)./H;    
        Output_Samples_block = ifft(ifftshift(Output_Freq_Samples_block));
        
        Output_Samples = [Output_Samples,Output_Samples_block(M:end)];
        
    end
end


function Output_Samples = CD_Compensation(Input_Samples, DTime, Fiber_Length, Beta2)
 
    Number_of_Samples = length(Input_Samples);
    MaxFreq = 0.5/DTime;
    DFreq = 2*MaxFreq/(Number_of_Samples-1);
    
    VFreq = (-1:Number_of_Samples-2)*DFreq;
    VFreq = VFreq-0.5*max(VFreq);
    VOmeg = 2*pi*VFreq;	
    
    %Dispersion Parameter
    Disper = (1j/2)*Beta2*VOmeg.^2*Fiber_Length; 
    % with - dispersion 
    % w/o - disp compensation
    % fft exp(-jwt) phase: jw*t-jbeta*z, engineer convention
  
    %dispersion
    Freq_Samples = fftshift(fft(Input_Samples.').');
    Output_Freq_Samples = Freq_Samples.*exp(Disper);
    Output_Samples = ifft(ifftshift(Output_Freq_Samples).').';
    
end


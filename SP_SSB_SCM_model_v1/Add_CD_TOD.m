function Output_Samples = Add_CD_TOD(Input_Samples, DTime, Fiber_Length, Beta1,Beta2,Beta3)
 
    Number_of_Samples = length(Input_Samples);
    MaxFreq = 0.5/DTime;
    VFreq = linspace(-MaxFreq,MaxFreq,Number_of_Samples);
    VOmeg = 2*pi*VFreq;	
    
    %Dispersion Parameter
    Disper = -1j*Beta1*VOmeg*Fiber_Length...
             -(1j/2)*Beta2*VOmeg.^2*Fiber_Length...
             -(1j/6)*Beta3*VOmeg.^3*Fiber_Length; 
    % with - dispersion 
    % w/o - disp compensation
    % fft exp(-jwt) phase: jw*t-jbeta*z, engineer convention
  
    %dispersion
    Freq_Samples = fftshift(fft(Input_Samples.').');
    Output_Freq_Samples = Freq_Samples.*exp(Disper);
    Output_Samples = ifft(ifftshift(Output_Freq_Samples).').';
    
end


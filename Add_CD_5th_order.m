function Output_Samples = Add_CD_5th_order(Input_Samples, DTime, Fiber_Length, Beta1,Beta2,Beta3,Beta4,Beta5)
 
    Number_of_Samples = length(Input_Samples);
    MaxFreq = 0.5/DTime;
    
    
    VFreq = linspace(-MaxFreq,MaxFreq,Number_of_Samples);
    
    VOmeg = 2*pi*VFreq;	
    
    %Dispersion Parameter
    Disper = -1*1j*Beta1*VOmeg*Fiber_Length...
             +1j*1/2*(1j^2)*Beta2*VOmeg.^2*Fiber_Length...
             +1*1/6*(1j^3)*Beta3*VOmeg.^3*Fiber_Length...
             -1j*1/24*(1j^4)*Beta4*VOmeg.^4*Fiber_Length...
             -1*1/120*(1j^5)*Beta5*VOmeg.^5*Fiber_Length; 
    % with - dispersion 
    % w/o - disp compensation
    % fft exp(-jwt) phase: jw*t-jbeta*z, engineer convention
  
    %dispersion
    Freq_Samples = fftshift(fft(Input_Samples.').');
    Output_Freq_Samples = Freq_Samples.*exp(Disper);
    Output_Samples = ifft(ifftshift(Output_Freq_Samples).').';
    
end


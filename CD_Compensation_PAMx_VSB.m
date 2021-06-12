% for VSB-PAMx only
function Output_Samples = CD_Compensation_PAMx_VSB(Input_Samples, DTime, ParamVSB,ParamFib,ParamSig)
    BaudRate = ParamSig.Baud_Rate;
    Fiber_Length = ParamFib.FiberLength;
    Beta2 = ParamFib.Beta2;
    Number_of_Samples = length(Input_Samples);
    MaxFreq = 0.5/DTime;
    DFreq = 2*MaxFreq/(Number_of_Samples-1);

    VFreq = (-1:Number_of_Samples-2)*DFreq;
    VFreq = VFreq-0.5*max(VFreq);
    VOmeg = 2*pi*VFreq;	

    %Dispersion Parameter
    Disper = -(1j/2)*Beta2*VOmeg.^2*Fiber_Length; 
    % with - dispersion 
    % w/o - disp compensation
    % fft exp(-jwt) phase: jw*t-jbeta*z, engineer convention

    %dispersion
    S1_F = fftshift(fft(Input_Samples));


    Freq_min = 0;
    Freq_max = BaudRate;


    
    Attenuation_dB = zeros(size(VFreq)); 

    Slope_Stop = ParamVSB.Opt_Flt_offset; 
    Slope_Start =Slope_Stop-ParamVSB.Opt_Flt_Suppression_dB/ParamVSB.Opt_Flt_Slope_dBper10GHz*10e9;
    if min(VFreq) < Slope_Start
        Attenuation_dB(VFreq<=Slope_Start)= -ParamVSB.Opt_Flt_Suppression_dB;
        length_slope = length(Attenuation_dB(VFreq>Slope_Start & VFreq <=Slope_Stop));
        Attenuation_dB(VFreq>Slope_Start & VFreq <=Slope_Stop) = linspace(-ParamVSB.Opt_Flt_Suppression_dB,0,length_slope);
    else
        ParamVSB.Opt_Flt_Suppression_dB = -min(VFreq)/10e9*ParamVSB.Opt_Flt_Slope_dBper10GHz;
        length_slope = length(Attenuation_dB(VFreq <=Slope_Stop));
        Attenuation_dB(VFreq <=Slope_Stop) = linspace(-ParamVSB.Opt_Flt_Suppression_dB,0,length_slope);
    end

    attenuation = 10.^(0.05*Attenuation_dB);
%     attenuation = fliplr(attenuation);
    
    H2 = fliplr(attenuation);


    H2 = H2.*exp(-1*Disper);
    H2(VFreq<=0) = 1;

    H1 = attenuation;
    H1 = H1.*exp(1*Disper);
    H1(VFreq<=0) = 1;



    S1_F = S1_F./(H2+H1);
    S1_F = S1_F.*(VFreq>Freq_min &VFreq <=Freq_max);


    Output_Samples = ifft(ifftshift(S1_F));

end
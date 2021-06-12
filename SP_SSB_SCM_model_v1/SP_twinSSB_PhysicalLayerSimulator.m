%% Tool Func
FT = inline('DTime*fftshift(fft(ifftshift(AtIn)))','AtIn','DTime'); %Analog Fourier transform

if ParamControl.new_capture
    %% DAC model (including Amp)
    % *********************************************************************
    if ParamControl.new_load
        [Tx_Time_Data_out1,Tx_Sequence1,ParamPS1] = zp_SP_SCM_SSB_Tx(ParamControl,ParamSig,ParamDAC,ParamFib);
        [Tx_Time_Data_out2,Tx_Sequence2,ParamPS2] = zp_SP_SCM_SSB_Tx(ParamControl,ParamSig,ParamDAC,ParamFib);
        Tx_Time_Data_out = Tx_Time_Data_out1 + conj(Tx_Time_Data_out2); 
    end
    % make longer Tx sequence
    Tx_Time_Data = [Tx_Time_Data_out,Tx_Time_Data_out];
    % Time-Freq Vectors
    DTime = 1/ParamDAC.DAC_Rate;
    Number_of_Samples = length(Tx_Time_Data);
    MaxFreq = 0.5/DTime; 
    
    DFreq = 2*MaxFreq/(Number_of_Samples-1);
    Freq_Vector = linspace(-MaxFreq,MaxFreq,Number_of_Samples);
    
    %% twin SSB always IQ
    if ParamControl.Preemphasis_or_Not
        switch ParamControl.DAC_LPF_option
            case 1
                Tx_Time_Data_I = real(anti_lpfilt(real(Tx_Time_Data),ParamDAC.DAC_Rate/1e9,...
                    ParamDAC.DAC_LPF_BW/1e9,ParamDAC.DAC_LPF_type,ParamDAC.DAC_LPF_order));
                Tx_Time_Data_Q = real(anti_lpfilt(imag(Tx_Time_Data),ParamDAC.DAC_Rate/1e9,...
                    ParamDAC.DAC_LPF_BW/1e9,ParamDAC.DAC_LPF_type,ParamDAC.DAC_LPF_order));
            case 2
                Tx_Time_Data = resample(Tx_Time_Data,88e9,ParamDAC.DAC_Rate);
                Tx_Time_Data_I = conv(real(Tx_Time_Data),ParamDAC.Preemphasis_FIR,'same');
                Tx_Time_Data_I = resample(Tx_Time_Data_I,ParamDAC.DAC_Rate,88e9);
                Tx_Time_Data_Q = conv(imag(Tx_Time_Data),ParamDAC.Preemphasis_FIR,'same');
                Tx_Time_Data_Q = resample(Tx_Time_Data_Q,ParamDAC.DAC_Rate,88e9);

        end
        Tx_Time_Data = Tx_Time_Data_I + 1i*Tx_Time_Data_Q;
        Tx_Time_Data = pwr_normalization(Tx_Time_Data);

        if ParamControl.Plot_Spectrum_or_Not
            figure;

            plot(Freq_Vector/1e9,max(10*log10(abs(FT(Tx_Time_Data(1:length(Freq_Vector)),DTime))),-150),'LineWidth',2);
            xlabel('Freq/GHz'); ylabel('Power/dBm');
            title('DAC waveform after Preemphasis');
        end

    end
    if ParamControl.Clipping_or_Not
        Tx_Time_Data_I = clipFunctionV2(real(Tx_Time_Data),ParamDAC.clipping_Prob,round(2^ParamDAC.qnbit_DAC));
        Tx_Time_Data_Q = clipFunctionV2(imag(Tx_Time_Data),ParamDAC.clipping_Prob,round(2^ParamDAC.qnbit_DAC));
        Tx_Time_Data = pwr_normalization(Tx_Time_Data_I+1i*Tx_Time_Data_Q);
    end
    % Add Quantization noise
    Tx_Time_Data_I = real(Tx_Time_Data);
    Tx_Time_Data_Q = imag(Tx_Time_Data);

    Tx_Time_Data_I = Quantization(Tx_Time_Data_I,ParamDAC.qnbit_DAC);
    Tx_Time_Data_Q = Quantization(Tx_Time_Data_Q,ParamDAC.qnbit_DAC);

    if ParamControl.DAC_LPF_or_Not
        switch ParamControl.DAC_LPF_option
            case 1
                Tx_Time_Data_I = real(lpfilt(Tx_Time_Data_I,ParamDAC.DAC_Rate/1e9,...
                    ParamDAC.DAC_LPF_BW/1e9,ParamDAC.DAC_LPF_type,ParamDAC.DAC_LPF_order));
                Tx_Time_Data_Q = real(lpfilt(Tx_Time_Data_Q,ParamDAC.DAC_Rate/1e9,...
                    ParamDAC.DAC_LPF_BW/1e9,ParamDAC.DAC_LPF_type,ParamDAC.DAC_LPF_order));
            case 2
                Tx_Time_Data_I = resample(Tx_Time_Data_I,88e9,ParamDAC.DAC_Rate);
                Tx_Time_Data_Q = resample(Tx_Time_Data_Q,88e9,ParamDAC.DAC_Rate);
                Tx_Time_Data_I = conv(Tx_Time_Data_I,ParamDAC.DAC_response_FIR,'same');
                Tx_Time_Data_I = resample(Tx_Time_Data_I,ParamDAC.DAC_Rate,88e9);
                Tx_Time_Data_Q = conv(Tx_Time_Data_Q,ParamDAC.DAC_response_FIR,'same');
                Tx_Time_Data_Q = resample(Tx_Time_Data_Q,ParamDAC.DAC_Rate,88e9);

        end
        Tx_Time_Data = Tx_Time_Data_I + 1i*Tx_Time_Data_Q;
        Tx_Time_Data = pwr_normalization(Tx_Time_Data);

        if ParamControl.Plot_Spectrum_or_Not
            figure;
            plot(Freq_Vector/1e9,max(10*log10(abs(FT(Tx_Time_Data(1:length(Freq_Vector)),DTime))),-150),'LineWidth',2);
            xlabel('Freq/GHz'); ylabel('Power/dBm');
            title('DAC waveform after DAC LPF');
        end
    else
        Tx_Time_Data = Tx_Time_Data_I + 1i*Tx_Time_Data_Q;
        Tx_Time_Data = pwr_normalization(Tx_Time_Data);
    end

    VI = pwr_normalization(Tx_Time_Data_I);
    VQ = pwr_normalization(Tx_Time_Data_Q);
    VM_I = max(max(VI),-min(VI));
    VM_Q = max(max(VQ),-min(VQ));

    ParamPhysicalModel.Vrms_over_VM_I = 1/VM_I;
    ParamPhysicalModel.Vrms_over_VM_Q = 1/VM_Q;
    Tx_Time_Data = pwr_normalization(Tx_Time_Data_I + 1i*Tx_Time_Data_Q);

    ParamPhysicalModel.Vm_pwr_normalized_I = max(max(real(Tx_Time_Data)),-min(real(Tx_Time_Data)));
    ParamPhysicalModel.Vm_pwr_normalized_Q = max(max(imag(Tx_Time_Data)),-min(imag(Tx_Time_Data)));

    
    if ParamControl.Add_DAC_noise
        IFFT_bin_length = ParamDAC.DAC_Rate/ParamSig.Baud_Rate;

        Tx_Time_Data_I = Add_AWGN_real(real(Tx_Time_Data),ParamDAC.DAC_SNR,IFFT_bin_length,1); 
        Tx_Time_Data_Q = Add_AWGN_real(imag(Tx_Time_Data),ParamDAC.DAC_SNR,IFFT_bin_length,1); 
        Tx_Time_Data = Tx_Time_Data_I + 1i*Tx_Time_Data_Q;

    end
    
    %% Modulator model
    % *********************************************************************
    % MZM nonlinearity
    if ParamControl.Modulator_Nonlinearity_or_Not && ParamControl.CSPR_tuning_case ~= 1       
        ParamPhysicalModel.Vm_over_Vpi = ParamSys.Vpp_over_Vpi/2;

        Tx_Time_Data_I = sin(pi/2*(real(Tx_Time_Data)/ParamPhysicalModel.Vm_pwr_normalized_I*ParamPhysicalModel.Vm_over_Vpi+...
                        ParamSys.Vbias_over_Vpi));
        Tx_Time_Data_Q = sin(pi/2*(imag(Tx_Time_Data)/ParamPhysicalModel.Vm_pwr_normalized_Q*ParamPhysicalModel.Vm_over_Vpi+...
                        ParamSys.Vbias_over_Vpi));
        Tx_Time_Data = Tx_Time_Data_I + 1i*Tx_Time_Data_Q;


        if ParamControl.CSPR_tuning_case == 2
            Modulator_path_pwr_ratio = 1 - ParamSys.Carrier_path_pwr_ratio; 
            Carrier = sqrt(ParamSys.Carrier_path_pwr_ratio)*sqrt(1/2);
            Tx_Time_Data = sqrt(Modulator_path_pwr_ratio)*sqrt(1/2)*Tx_Time_Data + Carrier;
        end
        Number_of_Samples = length(Tx_Time_Data);
        Power_Tx_Time_Data = norm(Tx_Time_Data)^2/Number_of_Samples;
        % Now the power of a full swing BPSK NRZ signal onto the modulator is 1
        ParamPhysicalModel.Tx_optical_Loss_dB = ParamMod.Modulator_Loss_dB - 10*log10(Power_Tx_Time_Data); 
    end
    % Modulator LPF
    if ParamControl.Modulator_LPF_or_Not
        Tx_Time_Data = lpfilt(Tx_Time_Data,ParamDAC.DAC_Rate,ParamMod.Modulator_LPF_BW,'gaussian',ParamMod.Modulator_LPF_order);
    end
    

    % add AWGN (laser + penalty)
    if ParamControl.Add_AWGN_beforeFiber_or_Not
        IFFT_bin_length = ParamDAC.DAC_Rate/ParamSig.Baud_Rate/2;
        Tx_Time_Data = Add_AWGN(Tx_Time_Data,ParamChan.SNR,IFFT_bin_length,1); 
    end
    

    
    %% laser model
    % *********************************************************************  
    % phase noise
    % add optical tone based on CSPR in the linear model
    if ParamControl.CSPR_tuning_case  == 1 || ParamControl.Curve_Fitting_or_Not
        Number_of_Samples = length(Tx_Time_Data);
        Power_Tx_Time_Data = norm(Tx_Time_Data)^2/Number_of_Samples;

        CSPR = 10^(ParamSys.CSPR_dB/10);
        Carrier = sqrt(Power_Tx_Time_Data*CSPR).*ones(size(Tx_Time_Data));

 
        if ParamControl.Add_PhaseNoise_or_Not
            Tx_Time_Data_resample = resample(Tx_Time_Data,1/1e-12,ParamDAC.DAC_Rate);
            Carrier_resample = resample(Carrier,1/1e-12,ParamDAC.DAC_Rate);
            Phase_signal = Add_PhaseNoise(size(Tx_Time_Data_resample),1/1e-12,ParamLas.Laser_Linewidth);
            Phase_carrier = circshift(Phase_signal,ParamSys.Path_Mismatch/1e-12,2);
            Tx_Time_Data_resample = Tx_Time_Data_resample.*exp(1i*Phase_signal);
            Carrier_resample = Carrier_resample.*exp(1i*Phase_carrier);
            Tx_Time_Data = resample(Tx_Time_Data_resample,ParamDAC.DAC_Rate,1/1e-12);
            Carrier = resample(Carrier_resample,ParamDAC.DAC_Rate,1/1e-12);
        end
        Tx_Time_Data = Tx_Time_Data + Carrier;     
    end
    
    if ParamControl.Add_PhaseNoise_or_Not && ParamControl.CSPR_tuning_case ~= 1
        Phase_signal = Add_PhaseNoise(size(Tx_Time_Data),ParamDAC.DAC_Rate,ParamLas.Laser_Linewidth);
        Tx_Time_Data = Tx_Time_Data.*exp(1i*Phase_signal);
    end
    
    if ParamControl.Plot_Spectrum_or_Not
        figure;
        plot(Freq_Vector/1e9,10*log10(abs(FT(Tx_Time_Data(1:length(Freq_Vector)),DTime))),'LineWidth',2);
        xlabel('Freq/GHz'); ylabel('Power/dBm');
        title('Modulated waveform');
    end    
    % Laser power
    if ParamControl.CSPR_tuning_case == 1 || ParamControl.Curve_Fitting_or_Not
        ParamPhysicalModel.LchPwr_dBm = ParamSys.Target_RxPwr_dBm + ParamChan.Fiber_attenuation + ParamChan.excess_loss_dB;
    else
        ParamPhysicalModel.LchPwr_dBm = ParamLas.laser_power_dBm -ParamPhysicalModel.Tx_optical_Loss_dB;
    end
    
    
    
    %% Generate Analog waveform
    % *********************************************************************
    ParamPhysicalModel.P0_dBm = ParamPhysicalModel.LchPwr_dBm - ParamChan.excess_loss_dB/2; 
    ParamPhysicalModel.P0 = 10^(0.1*ParamPhysicalModel.P0_dBm)*1e-3; 
    
    Tx_Time_Data = resample(Tx_Time_Data,ParamSys.Analog_Rate,ParamDAC.DAC_Rate);
    Tx_Time_Data = pwr_normalization(Tx_Time_Data)*sqrt(ParamPhysicalModel.P0);
    
    ParamPhysicalModel.Received_Pwr_dBm = ParamPhysicalModel.LchPwr_dBm - ParamChan.Fiber_attenuation - ParamChan.excess_loss_dB; 
    if ParamControl.Curve_Fitting_or_Not ||  ParamControl.CSPR_tuning_case == 1 ...
        || ParamControl.RxPreAmp_or_Not
        ParamPhysicalModel.Received_Pwr_dBm  = ParamSys.Target_RxPwr_dBm;
    end
    %% split step Fourier Transform (SSFT)    
    % *********************************************************************
    ParamPhysicalModel.AccumulateDistance = 0;
    dZ = ParamFib.FiberLength/ParamFib.numsec_SSFT/2;
    for idx = 1:ParamFib.numsec_SSFT 
        if ParamControl.Add_NL_or_Not
            ParamFib.Leff = (1-exp(-ParamFib.Loss_alpha*dZ))/ParamFib.Loss_alpha;
            NL_phase = -1i*ParamFib.gamma*ParamFib.Leff*...
                       abs(Tx_Time_Data).^2;
            Tx_Time_Data  = Tx_Time_Data .*exp(NL_phase).*exp(-ParamFib.Loss_alpha/2*dZ);
        end

        if ParamControl.Add_CD_or_Not
            DTime = 1/ParamSys.Analog_Rate;
            switch ParamControl.CD_order
                case 2
                    Tx_Time_Data = Add_CD(Tx_Time_Data,...
                            DTime, dZ*2, ...
                            ParamFib.Beta2_ref);
                case 3
                    Tx_Time_Data = Add_CD_TOD(Tx_Time_Data,...
                            DTime, dZ*2, ParamFib.Beta1_ref,...
                            ParamFib.Beta2_ref,ParamFib.Beta3_ref);
                case 5
                    Tx_Time_Data = Add_CD_5th_order(Tx_Time_Data,...
                            DTime, dZ*2, ParamFib.Beta1_ref,...
                            ParamFib.Beta2_ref,ParamFib.Beta3_ref,...
                            ParamFib.Beta4_ref,ParamFib.Beta5_ref);
            end
        end

        if ParamControl.Add_NL_or_Not
            ParamFib.Leff = (1-exp(-ParamFib.Loss_alpha*dZ))/ParamFib.Loss_alpha;
            NL_phase = -1i*ParamFib.gamma*ParamFib.Leff*...
                       abs(Tx_Time_Data ).^2;
            Tx_Time_Data  = Tx_Time_Data .*exp(NL_phase).*exp(-ParamFib.Loss_alpha/2*dZ);
        end
    end
    
    if ParamControl.Add_AWGN_afterFiber_or_Not
        
        IFFT_bin_length = ParamSys.Analog_Rate/ParamSig.Baud_Rate/2;
        Tx_Time_Data = Add_AWGN(Tx_Time_Data,ParamChan.SNR,IFFT_bin_length,1); 
    end
end
% Time-Freq Vectors update
DTime = 1/ParamSys.Analog_Rate;
Number_of_Samples = length(Tx_Time_Data);
MaxFreq = 0.5/DTime; 

DFreq = 2*MaxFreq/(Number_of_Samples-1);
Freq_Vector = linspace(-MaxFreq,MaxFreq,Number_of_Samples);
% at the receiver
%% optical filter

Tx_Freq_Data = fftshift(fft(Tx_Time_Data));
switch ParamControl.OBPF_option
    case 0
        Opt_Flt_Suppression = 10^(ParamVSB.Opt_Flt_Suppression_dB/10);
        Tx_Freq_Data(round(length(Tx_Freq_Data)/2):end) = ...
            Tx_Freq_Data(round(length(Tx_Freq_Data)/2):end).*Opt_Flt_Suppression;
    case 1
        Attenuation_dB1 = zeros(size(Freq_Vector)); 
        Slope_Stop1 = ParamVSB.Opt_Flt_offset1; 
        Slope_Start1 =Slope_Stop1-ParamVSB.Opt_Flt_Suppression_dB1/ParamVSB.Opt_Flt_Slope_dBper10GHz1*10e9;
        if min(Freq_Vector) < Slope_Start1
            Attenuation_dB1(Freq_Vector<=Slope_Start1)= -ParamVSB.Opt_Flt_Suppression_dB1;
            length_slope1 = length(Attenuation_dB1(Freq_Vector>Slope_Start1 & Freq_Vector <=Slope_Stop1));
            Attenuation_dB1(Freq_Vector>Slope_Start1 & Freq_Vector <=Slope_Stop1) = linspace(-ParamVSB.Opt_Flt_Suppression_dB1,0,length_slope1);
        else
            ParamVSB.Opt_Flt_Suppression_dB1 = (Slope_Stop1-min(Freq_Vector))/10e9*ParamVSB.Opt_Flt_Slope_dBper10GHz1;
            length_slope1 = length(Attenuation_dB1(Freq_Vector <=Slope_Stop1));
            Attenuation_dB1(Freq_Vector <=Slope_Stop1) = linspace(-ParamVSB.Opt_Flt_Suppression_dB1,0,length_slope1);
        end
        % Attenuation_dB = flipud(Attenuation_dB);
        if ParamControl.Plot_opt_filter_or_Not 
            figure; plot(Freq_Vector,Attenuation_dB1);
        end
        Attenuation1 = 10.^(0.05*Attenuation_dB1);
        if size(Attenuation1,1)>1
            flipAttenuation1 = flipud(Attenuation1);
        else
            flipAttenuation1 = fliplr(Attenuation1);
        end
        
        Attenuation_dB2 = zeros(size(Freq_Vector)); 
        Slope_Stop2 = ParamVSB.Opt_Flt_offset2; 
        Slope_Start2 =Slope_Stop2-ParamVSB.Opt_Flt_Suppression_dB2/ParamVSB.Opt_Flt_Slope_dBper10GHz2*10e9;
        if min(Freq_Vector) < Slope_Start2
            Attenuation_dB2(Freq_Vector<=Slope_Start2)= -ParamVSB.Opt_Flt_Suppression_dB2;
            length_slope2 = length(Attenuation_dB2(Freq_Vector>Slope_Start2 & Freq_Vector <=Slope_Stop2));
            Attenuation_dB2(Freq_Vector>Slope_Start2 & Freq_Vector <=Slope_Stop2) = linspace(-ParamVSB.Opt_Flt_Suppression_dB2,0,length_slope2);
        else
            ParamVSB.Opt_Flt_Suppression_dB2 = (Slope_Stop2-min(Freq_Vector))/10e9*ParamVSB.Opt_Flt_Slope_dBper10GHz2;
            length_slope2 = length(Attenuation_dB2(Freq_Vector <=Slope_Stop2));
            Attenuation_dB2(Freq_Vector <=Slope_Stop2) = linspace(-ParamVSB.Opt_Flt_Suppression_dB2,0,length_slope2);
        end
        % Attenuation_dB = flipud(Attenuation_dB);
        if ParamControl.Plot_opt_filter_or_Not 
            figure; plot(Freq_Vector,Attenuation_dB2);
        end
        Attenuation2 = 10.^(0.05*Attenuation_dB2);
        if size(Attenuation2,1)>1
            flipAttenuation2 = flipud(Attenuation2);
        else
            flipAttenuation2 = fliplr(Attenuation2);
        end
        
        Tx_Freq_Data1 = Tx_Freq_Data.*Attenuation1;
        Tx_Freq_Data2 = Tx_Freq_Data.*flipAttenuation2;
        
    case 2
        Freq_shift = 6e9;
        Tx_Freq_Data = Practical_OBPF(Tx_Freq_Data,Freq_Vector,Freq_shift);
    case 3
        Freq_shift = -5e9; % 1- -5e9 2- 0e9
        flt_opt = 1; % 1- 80 period 2- 112 period
        Tx_Freq_Data = SiPh_OBPF(Tx_Freq_Data,Freq_Vector,Freq_shift,flt_opt);
end
Tx_Time_Data1 = ifft(ifftshift(Tx_Freq_Data1));  
Tx_Time_Data2 = ifft(ifftshift(Tx_Freq_Data2));  

if ParamControl.Measure_CSPR_or_Not
    Tx_Freq_Data = fftshift(fft(Tx_Time_Data1));
    Carrier_Freq_Data = Tx_Freq_Data ; 
    BW_narrow_filter = 0.1e9;
    Carrier_Freq_Data(Freq_Vector < -BW_narrow_filter) = 0;
    Carrier_Freq_Data(Freq_Vector > BW_narrow_filter) = 0;
    Carrier_Time_Data = ifft(ifftshift(Carrier_Freq_Data));
    Number_of_Samples = length(Tx_Time_Data);
    Power_Tx_Time_Data = norm(Tx_Time_Data)^2/Number_of_Samples;
    Power_Carrier_Time_Data = norm(Carrier_Time_Data)^2/Number_of_Samples;
    Power_Signal_Time_Data = Power_Tx_Time_Data - Power_Carrier_Time_Data;
    ParamPhysicalModel.Measured_CSPR = 10*log10(Power_Carrier_Time_Data/Power_Signal_Time_Data);
end
if ParamControl.Plot_Spectrum_or_Not
    figure;
    plot(Freq_Vector/1e9,10*log10(abs(FT(Tx_Time_Data1,DTime))),'LineWidth',2);
    xlabel('Freq/GHz'); ylabel('Power/dBm');
    title('filtered waveform 1');
    figure;
    plot(Freq_Vector/1e9,10*log10(abs(FT(Tx_Time_Data2,DTime))),'LineWidth',2);
    xlabel('Freq/GHz'); ylabel('Power/dBm');
    title('filtered waveform 2');
end

    %% PD model
    % *********************************************************************
    Tx_Time_Data1 = pwr_normalization(Tx_Time_Data1);
    Tx_Time_Data2 = pwr_normalization(Tx_Time_Data2);
    Rx_I1 = Tx_Time_Data1.*conj(Tx_Time_Data1);
    Rx_I2 = Tx_Time_Data2.*conj(Tx_Time_Data2);
    % PD noise
    Received_Pwr = 10^(0.1*ParamPhysicalModel.Received_Pwr_dBm)*1e-3;%[W]
    if ParamControl.Add_PD_noise
        switch ParamControl.PD_case    
            case 1 % PIN
                delta_f = ParamPD.BW_pin;
                SNR_PD_Linear = (ParamPD.R_pin*Received_Pwr)^2/...
                                (2*ParamGen.q_electron*(ParamPD.R_pin*Received_Pwr+ParamPD.Id_pin)*delta_f+...
                                 4*ParamGen.kB*ParamGen.Temperature/ParamPD.RL_pin*ParamPD.Fn_pin*delta_f);
                Rx_I1 = Rx_I1*Received_Pwr*ParamPD.R_pin; 
                Rx_I2 = Rx_I2*Received_Pwr*ParamPD.R_pin; 
                
                ParamPhysicalModel.SNR_PD_dB = 10*log10(SNR_PD_Linear);
                Rx_I = resample(Rx_I,delta_f*2,ParamSys.Analog_Rate);
                IFFT_bin_length = 2*delta_f/ParamSig.Baud_Rate;
                Rx_I = Add_AWGN_real(Rx_I,ParamPhysicalModel.SNR_PD_dB,IFFT_bin_length,1);
                Rx_I = resample(Rx_I,ParamSys.Analog_Rate,delta_f*2);
                
                Rx_V = Rx_I*ParamPD.RL_pin;
                
            case 2 % APD
                delta_f = ParamPD.BW_apd;
                ParamPD.FA_apd = ParamPD.kA_apd*ParamPD.M_apd+(1-ParamPD.kA_apd)*(2-1/ParamPD.M_apd);
                SNR_PD_Linear = (ParamPD.M_apd*ParamPD.R_apd*Received_Pwr)^2/...
                                (2*ParamGen.q_electron*ParamPD.M_apd^2*ParamPD.FA_apd*(ParamPD.R_apd*Received_Pwr+ParamPD.Id_apd)*delta_f+...
                                 4*ParamGen.kB*ParamGen.Temperature/ParamPD.RL_apd*ParamPD.Fn_apd*delta_f);
                Rx_I = Rx_I*Received_Pwr*ParamPD.R_apd*ParamPD.M_apd; 
                
                ParamPhysicalModel.SNR_PD_dB = 10*log10(SNR_PD_Linear);
                Rx_I = resample(Rx_I,delta_f*2,ParamSys.Analog_Rate);
                IFFT_bin_length = 2*delta_f/ParamSig.Baud_Rate;
                Rx_I = Add_AWGN_real(Rx_I,ParamPhysicalModel.SNR_PD_dB,IFFT_bin_length,1);
                Rx_I = resample(Rx_I,ParamSys.Analog_Rate,delta_f*2);
                
                Rx_V = Rx_I*ParamPD.RL_apd;
                
                
            case 3 % PIN-TIA NEP model
                Rx_V1 = Rx_I1*Received_Pwr*ParamPD.CG; 
                Rx_V2 = Rx_I2*Received_Pwr*ParamPD.CG; 
                PD_noise_Vrms = ParamPD.CG*ParamPD.NEP; 
                ParamPhysicalModel.SNR_PD_dB = 20*log10(rms(Rx_V1)/PD_noise_Vrms);
                RL = 50;
                PD_noise = wgn(1,length(Rx_V1),10*log10(PD_noise_Vrms^2/RL),RL,'dBW','real');
                Rx_V1 = Rx_V1 + PD_noise;  
                Rx_V2 = Rx_V2 + PD_noise;  
        end
        ParamPD.SNR_PD_dB = ParamPhysicalModel.SNR_PD_dB;
    else
        Rx_V1 = Rx_I1*Received_Pwr*50;
        Rx_V2 = Rx_I2*Received_Pwr*50;
    end
    
    % PD freq response
    if ParamControl.PD_LPF_or_Not
        switch ParamControl.PD_case
            case 1
                PD_BW_ord = ParamPD.BW_PIN_order;
                PD_BW = ParamPD.BW_pin;
            case 2
                PD_BW_ord = ParamPD.BW_APD_order;
                PD_BW = ParamPD.BW_apd;
            case 3
                PD_BW_ord = ParamPD.BW_PIN_TIA_order;
                PD_BW = ParamPD.BW_PIN_TIA;
        end
        Rx_V1 = lpfilt(Rx_V1,ParamSys.Analog_Rate,PD_BW,'gaussian',PD_BW_ord);
        Rx_V2 = lpfilt(Rx_V2,ParamSys.Analog_Rate,PD_BW,'gaussian',PD_BW_ord);
    end
    
    %% ADC model
    % *********************************************************************
    % ADC noise (brick wall filter is also considered)
    if ParamControl.Add_ADC_thermal_noise
        Rx_V1 = resample(Rx_V1,ParamADC.ADC_BW*2,ParamSys.Analog_Rate);
        Rx_V2 = resample(Rx_V2,ParamADC.ADC_BW*2,ParamSys.Analog_Rate);
        ADC_noise_power_LN = 10^(ParamADC.ADC_noise_power_dBm/10)/12.5e9*ParamADC.ADC_BW;
        noise = wgn(1,length(Rx_V1),10*log10(ADC_noise_power_LN),50,'dBm','real');
        %ADC_noise_rms = rms(noise-mean(noise));
        %ADC_noise_Vpp = max(noise)-min(noise);
        Rx_V1 = Rx_V1 + noise;
        Rx_V2 = Rx_V2 + noise;
        ParamADC.SNR_ADC = 10*log10(get_pwr(Rx_V1)/get_pwr(noise));
        Rx_V1 = resample(Rx_V1,ParamSys.Analog_Rate,ParamADC.ADC_BW*2); 
        Rx_V2 = resample(Rx_V2,ParamSys.Analog_Rate,ParamADC.ADC_BW*2); 
    end
    
    % ADC detection
    Rx_V1 = resample(Rx_V1,ParamADC.ADC_Rate,ParamSys.Analog_Rate);
    Rx_V2 = resample(Rx_V2,ParamADC.ADC_Rate,ParamSys.Analog_Rate);
    
    ParamPhysicalModel.RV_rms_at_ADC = rms(Rx_V1-mean(Rx_V1));
    ParamPhysicalModel.RVpp_at_ADC = max(Rx_V1)-min(Rx_V1);
    
    % Rx_V = circshift(Rx_V,[0,2e3+2]);
    Rx_V1 = Quantization(real(Rx_V1),ParamADC.qnbit_ADC);
    Rx_V2 = Quantization(real(Rx_V2),ParamADC.qnbit_ADC);
    % Rx Side


    [ParamPhysicalModel.BER_avg, ParamPhysicalModel.BER_list, ParamPhysicalModel.SNR_list] = zp_SP_SCM_twinSSB_Rx(ParamControl,ParamRxDSP,ParamSig,ParamADC,ParamPS1,ParamPS2,...
                                            ParamFib,ParamMod,ParamPD,ParamVSB,...
                                            Rx_V1,Rx_V2,Tx_Sequence1,Tx_Sequence2);  



    
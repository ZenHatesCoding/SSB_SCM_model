%% Tool Func
FT = inline('DTime*fftshift(fft(ifftshift(AtIn)))','AtIn','DTime'); %Analog Fourier transform

if ParamControl.new_capture
    %% DAC model (including Amp)
    % *********************************************************************
    if ParamControl.new_load
        switch ParamControl.PAM_or_QAM
            case 1
                [Tx_Time_Data_out,Tx_Sequence,ParamPS] = zp_SP_SCM_DSB_Tx(ParamControl,ParamSig,ParamDAC,ParamFib);
                if ParamControl.VSB_or_Not == 0
                    Tx_Time_Data_out = hilbert(Tx_Time_Data_out);
                end
            case 2
                [Tx_Time_Data_out,Tx_Sequence,ParamPS] = zp_SP_SCM_SSB_Tx(ParamControl,ParamSig,ParamDAC,ParamFib);
        end
    end
    % make longer Tx sequence
%     Tx_Time_Data = [Tx_Time_Data_out,Tx_Time_Data_out];
    Tx_Time_Data = Tx_Time_Data_out;


    % Time-Freq Vectors
    DTime = 1/ParamDAC.DAC_Rate;
    Number_of_Samples = length(Tx_Time_Data);
    MaxFreq = 0.5/DTime; 
    
    DFreq = 2*MaxFreq/(Number_of_Samples-1);
    Freq_Vector = linspace(-MaxFreq,MaxFreq,Number_of_Samples);
    
    if ParamControl.VSB_or_Not      
        % add quantization noise
        if ParamControl.Quantization_or_Not
        Tx_Time_Data = Quantization(Tx_Time_Data,ParamDAC.qnbit_DAC);
        end
        if ParamControl.DAC_LPF_or_Not
            switch ParamControl.DAC_LPF_option
                case 1
                    Tx_Time_Data = real(lpfilt(Tx_Time_Data,ParamDAC.DAC_Rate/1e9,...
                                ParamDAC.DAC_LPF_BW/1e9,ParamDAC.DAC_LPF_type,ParamDAC.DAC_LPF_order));
                case 2
                    Tx_Time_Data = resample(Tx_Time_Data,88e9,ParamDAC.DAC_Rate);
                    Tx_Time_Data = conv(Tx_Time_Data,ParamDAC.DAC_response_FIR,'same');
                    Tx_Time_Data = resample(Tx_Time_Data,ParamDAC.DAC_Rate,88e9);
            end
            Tx_Time_Data = pwr_normalization(Tx_Time_Data);
            if ParamControl.Plot_Spectrum_or_Not 
                figure;
                plot(Freq_Vector/1e9,max(10*log10(abs(FT(Tx_Time_Data,DTime))),-150),'LineWidth',2);
                xlabel('Freq/GHz'); ylabel('Power/dBm');
                title('DAC waveform after DAC LPF');
            end
        end
        Tx_Time_Data = pwr_normalization(Tx_Time_Data);
        
        VM = max(max(Tx_Time_Data),-min(Tx_Time_Data));
        
        ParamPhysicalModel.Vrms_over_VM = 1/VM;
        ParamPhysicalModel.Vm_pwr_normalized = max(max(Tx_Time_Data),-min(Tx_Time_Data));
    else % IQ based SSB
        % Add Quantization noise
        Tx_Time_Data_I = real(Tx_Time_Data);
        Tx_Time_Data_Q = imag(Tx_Time_Data);
        if ParamControl.Quantization_or_Not
        Tx_Time_Data_I = Quantization(Tx_Time_Data_I,ParamDAC.qnbit_DAC);
        Tx_Time_Data_Q = Quantization(Tx_Time_Data_Q,ParamDAC.qnbit_DAC);
        end
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
                plot(Freq_Vector/1e9,max(20*log10(abs(FT(Tx_Time_Data(1:length(Freq_Vector)),DTime))),-150),'LineWidth',2);
                xlabel('Freq/GHz'); ylabel('Power/dBm');
                title('DAC waveform after DAC LPF');
            end
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
    end
    
    if ParamControl.Add_DAC_noise
        IFFT_bin_length = ParamDAC.DAC_Rate/ParamSig.Baud_Rate;
        if ParamControl.VSB_or_Not
            Tx_Time_Data = Add_AWGN_real(Tx_Time_Data,ParamDAC.DAC_SNR,IFFT_bin_length,1); 
        else
            Tx_Time_Data_I = Add_AWGN_real(real(Tx_Time_Data),ParamDAC.DAC_SNR,IFFT_bin_length,1); 
            Tx_Time_Data_Q = Add_AWGN_real(imag(Tx_Time_Data),ParamDAC.DAC_SNR,IFFT_bin_length,1); 
            Tx_Time_Data = Tx_Time_Data_I + 1i*Tx_Time_Data_Q;
        end
    end
    
  


    %% Modulator model
    % *********************************************************************
    % MZM nonlinearity
    if ParamControl.Modulator_Nonlinearity_or_Not && ParamControl.CSPR_tuning_case ~= 1       
        ParamPhysicalModel.Vm_over_Vpi = ParamSys.Vpp_over_Vpi/2;
        if ParamControl.VSB_or_Not
            Tx_Time_Data = sin(pi/2*(Tx_Time_Data/ParamPhysicalModel.Vm_pwr_normalized*ParamPhysicalModel.Vm_over_Vpi+...
                           ParamSys.Vbias_over_Vpi));
            % Modulator LPF
            if ParamControl.Modulator_LPF_or_Not
                Tx_Time_Data = lpfilt(Tx_Time_Data,ParamDAC.DAC_Rate,ParamMod.Modulator_LPF_BW,'gaussian',ParamMod.Modulator_LPF_order);
            end
            
        else
            switch ParamControl.MZM_option
                case 1
                    Tx_Time_Data_I = 1/2*sin(pi/2*(real(Tx_Time_Data)/ParamPhysicalModel.Vm_pwr_normalized_I*ParamPhysicalModel.Vm_over_Vpi+...
                                    ParamSys.Vbias_over_Vpi));
                    Tx_Time_Data_Q = 1/2*sin(pi/2*(imag(Tx_Time_Data)/ParamPhysicalModel.Vm_pwr_normalized_Q*ParamPhysicalModel.Vm_over_Vpi+...
                                    ParamSys.Vbias_over_Vpi));
                    
                    if ParamControl.Modulator_LPF_or_Not
                        Tx_Time_Data_I = lpfilt(Tx_Time_Data_I,ParamDAC.DAC_Rate,ParamMod.Modulator_LPF_BW,'gaussian',ParamMod.Modulator_LPF_order);
                        Tx_Time_Data_Q = lpfilt(Tx_Time_Data_Q,ParamDAC.DAC_Rate,ParamMod.Modulator_LPF_BW,'gaussian',ParamMod.Modulator_LPF_order);
                    end 
                    
                    IQ_Ampimbalance = 10^(-0.1*ParamMod.IQAmpImbalance_dB);
                    Tx_Time_Data_Q = Tx_Time_Data_Q*exp(1i*pi/180*ParamMod.IQPhaseImbalance_deg);
                    Tx_Time_Data_Q = Tx_Time_Data_Q*sqrt(2/(IQ_Ampimbalance^2+1));
                    omega = linspace(-ParamDAC.DAC_Rate/2,ParamDAC.DAC_Rate/2,length(Tx_Time_Data_Q));
                    Tx_Time_Data_Q = ifft(ifftshift(fftshift(fft(Tx_Time_Data_Q)).*exp(1i*omega*ParamMod.IQskew))); 
                    
                    Tx_Time_Data_I = Tx_Time_Data_I*sqrt((2*IQ_Ampimbalance^2)/(IQ_Ampimbalance^2+1));
                    
                    
                    Tx_Time_Data = Tx_Time_Data_I + 1i*Tx_Time_Data_Q;
                    ParamMod.Modulator_Loss_dB= ParamMod.Modulator_Loss_dB+2*ParamMod.Y_branch_Loss_dB;
                case 2
                    Tx_Time_Data = 1/2*(exp(1i*pi*(real(Tx_Time_Data)/ParamPhysicalModel.Vm_pwr_normalized_I*ParamPhysicalModel.Vm_over_Vpi+...
                                    ParamSys.Vbias_over_Vpi/2))+...
                                        exp(1i*pi*(imag(Tx_Time_Data)/ParamPhysicalModel.Vm_pwr_normalized_Q*ParamPhysicalModel.Vm_over_Vpi-...
                                    ParamSys.Vbias_over_Vpi/2)));
            end
        end

        if ParamControl.CSPR_tuning_case == 2
            Modulator_path_pwr_ratio = 1 - ParamSys.Carrier_path_pwr_ratio; 
            Carrier = sqrt(ParamSys.Carrier_path_pwr_ratio)*sqrt(1/2)*sqrt(10^(-0.1*ParamMod.carrier_path_loss_dB));
            
            

            Tx_Time_Data = sqrt(Modulator_path_pwr_ratio)*sqrt(1/2)*Tx_Time_Data*sqrt(10^(-0.1*ParamMod.Modulator_Loss_dB))...
                           + Carrier;
            Number_of_Samples = length(Tx_Time_Data);
            Power_Tx_Time_Data = sum(abs(Tx_Time_Data).^2)/Number_of_Samples;
            ParamPhysicalModel.Tx_optical_Loss_dB = -10*log10(Power_Tx_Time_Data)+2*ParamMod.Y_branch_Loss_dB+ParamMod.MMI_Loss_dB;
            
        else
            Number_of_Samples = length(Tx_Time_Data);
%             Power_Tx_Time_Data = norm(Tx_Time_Data)^2/Number_of_Samples;
            Power_Tx_Time_Data = sum(abs(Tx_Time_Data).^2)/Number_of_Samples;
            % Now the power of a full swing BPSK NRZ signal onto the modulator is 1
            ParamPhysicalModel.Tx_optical_Loss_dB = ParamMod.Modulator_Loss_dB - 10*log10(Power_Tx_Time_Data); 
        end
        ParamPhysicalModel.Tx_optical_Loss_dB = ParamPhysicalModel.Tx_optical_Loss_dB  + 2*ParamMod.SSC_Loss_dB;
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
    
    % Laser power
    if ParamControl.CSPR_tuning_case == 1 || ParamControl.Curve_Fitting_or_Not
        ParamPhysicalModel.LchPwr_dBm = ParamSys.Target_RxPwr_dBm + ParamChan.Fiber_attenuation + 3*ParamChan.coupler_loss_dB;
    else
        ParamPhysicalModel.LchPwr_dBm = ParamLas.laser_power_dBm -ParamPhysicalModel.Tx_optical_Loss_dB-2*ParamChan.coupler_loss_dB;
    end
    
    
    if ParamControl.Plot_Spectrum_or_Not
        figure;
        plot(Freq_Vector/1e9,20*log10(abs(FT(Tx_Time_Data(1:length(Freq_Vector)),DTime))),'LineWidth',2);
        xlabel('Freq/GHz'); ylabel('Power/dBm');
        title('Modulated waveform');
    end
    % add AWGN (laser + penalty)
    if ParamControl.Add_AWGN_beforeFiber_or_Not
        IFFT_bin_length = ParamDAC.DAC_Rate/ParamSig.Baud_Rate;
        Tx_Time_Data = Add_AWGN(Tx_Time_Data,ParamChan.SNR,IFFT_bin_length,1); 
    end


    
    if ParamControl.VSB_or_Not && ParamControl.VSB_at_Tx_or_Rx == 1
        ParamPhysicalModel.Opt_Flt_loss_dB = ParamVSB.Opt_Flt_ILoss_dB;
        Power_before_Opt_Flt = get_pwr(Tx_Time_Data);
        Tx_Freq_Data = fftshift(fft(Tx_Time_Data));
        switch ParamControl.OBPF_option
            case 0
                Opt_Flt_Suppression = 10^(ParamVSB.Opt_Flt_Suppression_dB/10);
                Tx_Freq_Data(round(length(Tx_Freq_Data)/2):end) = ...
                    Tx_Freq_Data(round(length(Tx_Freq_Data)/2):end).*Opt_Flt_Suppression;
            case 1
                Attenuation_dB = zeros(size(Freq_Vector)); 
                Slope_Stop = ParamVSB.Opt_Flt_offset+ParamVSB.Opt_Flt_drift; 
                Slope_Start =Slope_Stop-ParamVSB.Opt_Flt_Suppression_dB/ParamVSB.Opt_Flt_Slope_dBper10GHz*10e9;
                if min(Freq_Vector) < Slope_Start
                    Attenuation_dB(Freq_Vector<=Slope_Start)= -ParamVSB.Opt_Flt_Suppression_dB;
                    length_slope = length(Attenuation_dB(Freq_Vector>Slope_Start & Freq_Vector <=Slope_Stop));
                    Attenuation_dB(Freq_Vector>Slope_Start & Freq_Vector <=Slope_Stop) = linspace(-ParamVSB.Opt_Flt_Suppression_dB,0,length_slope);
                else
                    ParamVSB.Opt_Flt_Suppression_dB = (Slope_Stop-min(Freq_Vector))/10e9*ParamVSB.Opt_Flt_Slope_dBper10GHz;
                    length_slope = length(Attenuation_dB(Freq_Vector <=Slope_Stop));
                    Attenuation_dB(Freq_Vector <=Slope_Stop) = linspace(-ParamVSB.Opt_Flt_Suppression_dB,0,length_slope);
                end

                if ParamControl.Plot_VSB_filter_or_Not 
                    figure; plot(Freq_Vector,Attenuation_dB);
                end
                Attenuation = 10.^(0.05*Attenuation_dB);
                if ParamControl.LSB_or_RSB
                    Attenuation = fliplr(Attenuation);
                end
                Tx_Freq_Data = Tx_Freq_Data.*Attenuation;
            case 2
                Freq_shift = 6e9;
                Tx_Freq_Data = Practical_OBPF(Tx_Freq_Data,Freq_Vector,Freq_shift);
            case 3
                Freq_shift = -5e9; % 1- -5e9 2- 0e9
                flt_opt = 1; % 1- 80 period 2- 112 period
                Tx_Freq_Data = SiPh_OBPF(Tx_Freq_Data,Freq_Vector,Freq_shift,flt_opt);
            case 4
                Tx_Freq_Data(1:(length(Tx_Freq_Data)/2-1))=0;
            case 5
                Freq_shift = ParamVSB.Opt_Flt_offset+ParamVSB.Opt_Flt_drift;               
                Tx_Freq_Data = Lumentum_OBPF(Tx_Freq_Data,Freq_Vector,Freq_shift,ParamControl);
        end
        Tx_Time_Data = ifft(ifftshift(Tx_Freq_Data));  
        % loss of filter
        Power_after_Opt_Flt = get_pwr(Tx_Time_Data);
        ParamPhysicalModel.Opt_Flt_loss_dB = ParamPhysicalModel.Opt_Flt_loss_dB - 10*log10(Power_after_Opt_Flt/Power_before_Opt_Flt);
        ParamPhysicalModel.LchPwr_dBm  = ParamPhysicalModel.LchPwr_dBm -...
                                              ParamPhysicalModel.Opt_Flt_loss_dB;


    end
    


    
    
    %% Generate Analog waveform
    % *********************************************************************
    ParamPhysicalModel.P0_dBm = ParamPhysicalModel.LchPwr_dBm ; 
    ParamPhysicalModel.P0 = 10^(0.1*ParamPhysicalModel.P0_dBm)*1e-3; 
    
    Tx_Time_Data = resample(Tx_Time_Data,ParamSys.Analog_Rate,ParamDAC.DAC_Rate);
    Tx_Time_Data = pwr_normalization(Tx_Time_Data)*sqrt(ParamPhysicalModel.P0);
    
    ParamPhysicalModel.Received_Pwr_dBm = ParamPhysicalModel.LchPwr_dBm - ParamChan.Fiber_attenuation - 1*ParamChan.coupler_loss_dB-...
                                          ParamChan.extra_loss_dB; 
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
                            DTime, dZ*2, ParamFib.Beta1_ref-ParamFib.Beta1_ref,...
                            ParamFib.Beta2_ref,ParamFib.Beta3_ref);
                case 5
                    Tx_Time_Data = Add_CD_5th_order(Tx_Time_Data,...
                            DTime, dZ*2, ParamFib.Beta1_ref-ParamFib.Beta1_ref,...
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
        IFFT_bin_length = ParamSys.Analog_Rate/ParamSig.Baud_Rate;
        Tx_Time_Data = Add_AWGN(Tx_Time_Data,ParamChan.SNR,IFFT_bin_length,1); 
    end
    % update Freq Vector
    DTime = 1/ParamSys.Analog_Rate;
    Number_of_Samples = length(Tx_Time_Data);
    MaxFreq = 0.5/DTime; 

    DFreq = 2*MaxFreq/(Number_of_Samples-1);
    Freq_Vector = linspace(-MaxFreq,MaxFreq,Number_of_Samples);
%     Tx_Time_Data = pwr_normalization(Tx_Time_Data);

    if ParamControl.VSB_or_Not && ParamControl.VSB_at_Tx_or_Rx == 0
        ParamPhysicalModel.Opt_Flt_loss_dB = ParamVSB.Opt_Flt_ILoss_dB;
        Power_before_Opt_Flt = get_pwr(Tx_Time_Data);
        Tx_Freq_Data = fftshift(fft(Tx_Time_Data));
        switch ParamControl.OBPF_option
            case 0
                Opt_Flt_Suppression = 10^(ParamVSB.Opt_Flt_Suppression_dB/10);
                Tx_Freq_Data(round(length(Tx_Freq_Data)/2):end) = ...
                    Tx_Freq_Data(round(length(Tx_Freq_Data)/2):end).*Opt_Flt_Suppression;
            case 1
                Attenuation_dB = zeros(size(Freq_Vector)); 
                Slope_Stop = ParamVSB.Opt_Flt_offset+ParamVSB.Opt_Flt_drift; 
                Slope_Start =Slope_Stop-ParamVSB.Opt_Flt_Suppression_dB/ParamVSB.Opt_Flt_Slope_dBper10GHz*10e9;
                if min(Freq_Vector) < Slope_Start
                    Attenuation_dB(Freq_Vector<=Slope_Start)= -ParamVSB.Opt_Flt_Suppression_dB;
                    length_slope = length(Attenuation_dB(Freq_Vector>Slope_Start & Freq_Vector <=Slope_Stop));
                    Attenuation_dB(Freq_Vector>Slope_Start & Freq_Vector <=Slope_Stop) = linspace(-ParamVSB.Opt_Flt_Suppression_dB,0,length_slope);
                else
                    ParamVSB.Opt_Flt_Suppression_dB = (Slope_Stop-min(Freq_Vector))/10e9*ParamVSB.Opt_Flt_Slope_dBper10GHz;
                    length_slope = length(Attenuation_dB(Freq_Vector <=Slope_Stop));
                    Attenuation_dB(Freq_Vector <=Slope_Stop) = linspace(-ParamVSB.Opt_Flt_Suppression_dB,0,length_slope);
                end
                
                if ParamControl.Plot_VSB_filter_or_Not 
                    figure; plot(Freq_Vector,Attenuation_dB);
                end
                Attenuation = 10.^(0.05*Attenuation_dB);
                if ParamControl.LSB_or_RSB
                    Attenuation = fliplr(Attenuation);
                end
                Tx_Freq_Data = Tx_Freq_Data.*Attenuation;
            case 2
                Freq_shift = 6e9;
                Tx_Freq_Data = Practical_OBPF(Tx_Freq_Data,Freq_Vector,Freq_shift);
            case 3
                Freq_shift = -5e9; % 1- -5e9 2- 0e9
                flt_opt = 1; % 1- 80 period 2- 112 period
                Tx_Freq_Data = SiPh_OBPF(Tx_Freq_Data,Freq_Vector,Freq_shift,flt_opt);
            case 4
                Tx_Freq_Data(1:(length(Tx_Freq_Data)/2-1))=0;
            case 5
                Freq_shift = ParamVSB.Opt_Flt_offset+ParamVSB.Opt_Flt_drift;
                Tx_Freq_Data = Lumentum_OBPF(Tx_Freq_Data,Freq_Vector,Freq_shift,ParamControl);
        end
        Tx_Time_Data = ifft(ifftshift(Tx_Freq_Data));  
        
        % loss of filter
        Power_after_Opt_Flt = get_pwr(Tx_Time_Data);
        ParamPhysicalModel.Opt_Flt_loss_dB = ParamPhysicalModel.Opt_Flt_loss_dB - 10*log10(Power_after_Opt_Flt/Power_before_Opt_Flt);
        ParamPhysicalModel.Received_Pwr_dBm = ParamPhysicalModel.Received_Pwr_dBm-...
                                              ParamPhysicalModel.Opt_Flt_loss_dB;

    end
    if ParamControl.Plot_Spectrum_or_Not
        figure;
        
        plot(Freq_Vector/1e9,20*log10(abs(FT(Tx_Time_Data,DTime))),'LineWidth',2);
        xlabel('Freq/GHz'); ylabel('Power/dBm');
        title('received waveform');
    end       
    
    
end
   % measure CSPR
    if ParamControl.Measure_CSPR_or_Not
        Tx_Freq_Data = fftshift(fft(Tx_Time_Data));
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
    




    %% PD model
    % *********************************************************************
    Tx_Time_Data = pwr_normalization(Tx_Time_Data);
    Rx_I = Tx_Time_Data.*conj(Tx_Time_Data);
    % trick
%     Rx_I = Tx_Time_Data;
    
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
        Rx_I = lpfilt(Rx_I,ParamSys.Analog_Rate,PD_BW,'gaussian',PD_BW_ord);
    end
    
    % PD noise
    
    Received_Pwr = 10^(0.1*ParamPhysicalModel.Received_Pwr_dBm)*1e-3;%[W]
    if ParamControl.Add_PD_noise
        switch ParamControl.PD_case    
            case 1 % PIN
                delta_f = ParamPD.BW_pin;
                SNR_PD_Linear = (ParamPD.R_pin*Received_Pwr)^2/...
                                (2*ParamGen.q_electron*(ParamPD.R_pin*Received_Pwr+ParamPD.Id_pin)*delta_f+...
                                 4*ParamGen.kB*ParamGen.Temperature/ParamPD.RL_pin*ParamPD.Fn_pin*delta_f);
                Rx_I = Rx_I*Received_Pwr*ParamPD.R_pin; 
                
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
                Rx_I = resample(Rx_I,1,2); % trick
                Rx_V = Rx_I*Received_Pwr*ParamPD.CG; 
                PD_noise_Vrms = ParamPD.CG*ParamPD.NEP; 
                ParamPhysicalModel.SNR_PD_dB = 20*log10(rms(Rx_V)/PD_noise_Vrms);
                RL = 50;
                PD_noise = wgn(1,length(Rx_V),10*log10(PD_noise_Vrms^2/RL),RL,'dBW','real');
                Rx_V = Rx_V + PD_noise;               
        end
        ParamPD.SNR_PD_dB = ParamPhysicalModel.SNR_PD_dB;
    else
        Rx_V = Rx_I*Received_Pwr*50;
    end
    

    
    %% ADC model
    % *********************************************************************
    % ADC noise (brick wall filter is also considered)
    if ParamControl.Add_ADC_thermal_noise
        Rx_V = resample(Rx_V,ParamADC.ADC_BW*2,ParamSys.Analog_Rate/2);
        ADC_noise_power_LN = 10^(ParamADC.ADC_noise_power_dBm/10)/12.5e9*ParamADC.ADC_BW;
        noise = wgn(1,length(Rx_V),10*log10(ADC_noise_power_LN),50,'dBm','real');
        %ADC_noise_rms = rms(noise-mean(noise));
        %ADC_noise_Vpp = max(noise)-min(noise);
        Rx_V = Rx_V + noise;
        ParamADC.SNR_ADC = 10*log10(get_pwr(Rx_V)/get_pwr(noise));
        Rx_V = resample(Rx_V,ParamSys.Analog_Rate/2,ParamADC.ADC_BW*2); 
        ParamPhysicalModel.SNR_ADC =  ParamADC.SNR_ADC ;
    end
    
    % ADC detection
    Rx_V = resample(Rx_V,ParamADC.ADC_Rate,ParamSys.Analog_Rate/2);
    
    ParamPhysicalModel.RV_rms_at_ADC = rms(Rx_V-mean(Rx_V));
    ParamPhysicalModel.RVpp_at_ADC = max(Rx_V)-min(Rx_V);
    
    % Rx_V = circshift(Rx_V,[0,2e3+2]);
    if ParamControl.Quantization_or_Not
        Rx_V = Quantization(real(Rx_V),ParamADC.qnbit_ADC);
    end
    
    
    
    % Rx Side
    [ParamPhysicalModel.BER_avg, ParamPhysicalModel.BER_list, ParamPhysicalModel.SNR_list] = zp_SP_SCM_SSB_Rx(ParamControl,ParamRxDSP,ParamSig,ParamADC,ParamPS,ParamFib,ParamMod,ParamPD,ParamVSB,Rx_V,Tx_Sequence);  
    
    %% display on the screen
%     disp(ParamPhysicalModel);

    
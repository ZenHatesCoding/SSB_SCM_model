%% Tool Func
FT = inline('DTime*fftshift(fft(ifftshift(AtIn)))','AtIn','DTime'); %Analog Fourier transform

tic

if ParamControl.New_Simulation
    for numChannel = 1:ParamSig.numChannel
        %% DAC model (including Amp)
        % *********************************************************************

        [Tx_Time_Data_out,Tx_Sequence,ParamPS] = zp_SP_SCM_DSB_Tx(ParamControl,ParamSig,ParamDAC,ParamFib);
        Tx_Sequence_CWDM{numChannel} = Tx_Sequence;
        % make longer Tx sequence
        Tx_Time_Data = [Tx_Time_Data_out,Tx_Time_Data_out];
        % Time-Freq Vectors
        DTime = 1/ParamDAC.DAC_Rate;
        Number_of_Samples = length(Tx_Time_Data);
        MaxFreq = 0.5/DTime; 

        DFreq = 2*MaxFreq/(Number_of_Samples-1);
        Freq_Vector = linspace(-MaxFreq,MaxFreq,Number_of_Samples);

        if ParamControl.IM_or_Not      
            if ParamControl.Preemphasis_or_Not
                switch ParamControl.DAC_LPF_option
                    case 1
                        Tx_Time_Data = real(anti_lpfilt(Tx_Time_Data,ParamDAC.DAC_Rate/1e9,...
                            ParamDAC.DAC_LPF_BW/1e9,ParamDAC.DAC_LPF_type,ParamDAC.DAC_LPF_order));

                    case 2
                        Tx_Time_Data = resample(Tx_Time_Data,88e9,ParamDAC.DAC_Rate);
                        Tx_Time_Data = conv(Tx_Time_Data,ParamDAC.Preemphasis_FIR,'same');
                        Tx_Time_Data = resample(Tx_Time_Data,ParamDAC.DAC_Rate,88e9);
                end
                Tx_Time_Data = pwr_normalization(Tx_Time_Data);
                if ParamControl.Plot_Spectrum_or_Not 
                    figure;
                    plot(Freq_Vector/1e9,max(10*log10(abs(FT(Tx_Time_Data,DTime))),-150),'LineWidth',2);
                    xlabel('Freq/GHz'); ylabel('Power/dBm');
                    title('DAC waveform after Pre-emphasis');
                end
            end

            if ParamControl.Clipping_or_Not
                Tx_Time_Data = clipFunctionV2(Tx_Time_Data,ParamDAC.clipping_Prob,round(2^ParamDAC.qnbit_DAC));
                Tx_Time_Data = pwr_normalization(Tx_Time_Data);
            end
            % add quantization noise
            Tx_Time_Data = Quantization(Tx_Time_Data,ParamDAC.qnbit_DAC);

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
            if ParamControl.IM_or_Not
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
        if ParamControl.Modulator_Nonlinearity_or_Not     
            ParamPhysicalModel.Vm_over_Vpi = ParamSys.Vpp_over_Vpi/2;
            if ParamControl.IM_or_Not
                Tx_Time_Data = sin(pi/2*(Tx_Time_Data/ParamPhysicalModel.Vm_pwr_normalized*ParamPhysicalModel.Vm_over_Vpi+...
                               ParamSys.Vbias_over_Vpi));

            else
                Tx_Time_Data_I = sin(pi/2*(real(Tx_Time_Data)/ParamPhysicalModel.Vm_pwr_normalized_I*ParamPhysicalModel.Vm_over_Vpi+...
                                ParamSys.Vbias_over_Vpi));
                Tx_Time_Data_Q = sin(pi/2*(imag(Tx_Time_Data)/ParamPhysicalModel.Vm_pwr_normalized_Q*ParamPhysicalModel.Vm_over_Vpi+...
                                ParamSys.Vbias_over_Vpi));
                Tx_Time_Data = Tx_Time_Data_I + 1i*Tx_Time_Data_Q;

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

        if ParamControl.Plot_Spectrum_or_Not
            figure;
            plot(Freq_Vector/1e9,10*log10(abs(FT(Tx_Time_Data(1:length(Freq_Vector)),DTime))),'LineWidth',2);
            xlabel('Freq/GHz'); ylabel('Power/dBm');
            title('Modulated waveform');
        end

        %% laser model
        % *********************************************************************  
        % phase noise
        % add optical tone based on CSPR in the linear model
        if ParamControl.Add_PhaseNoise_or_Not 
            Phase_signal = Add_PhaseNoise(size(Tx_Time_Data),ParamDAC.DAC_Rate,ParamLas.Laser_Linewidth);
            Tx_Time_Data = Tx_Time_Data.*exp(1i*Phase_signal);
        end

        % Laser power

        ParamPhysicalModel.LchPwr_dBm = ParamLas.laser_power_dBm -ParamPhysicalModel.Tx_optical_Loss_dB;


        % add AWGN (laser + penalty)
        if ParamControl.Add_AWGN_or_Not
            IFFT_bin_length = ParamDAC.DAC_Rate/ParamSig.Baud_Rate;
            Tx_Time_Data = Add_AWGN(Tx_Time_Data,ParamChan.SNR,IFFT_bin_length,1); 
        end

        Tx_Time_Data_CWDM{numChannel} = pwr_normalization(Tx_Time_Data);
    end 
    

    % *************************************************************************
    % *************************************************************************
    % *************************************************************************
    if ParamControl.Coupled_NLS_or_Not == 0
        %% Generate Analog waveform
        % *********************************************************************
        Tx_Time_Data_CWDM_mux = 0;
        ParamPhysicalModel.P0_dBm = ParamPhysicalModel.LchPwr_dBm - ParamChan.excess_loss_dB/2; 
        ParamPhysicalModel.P0 = 10^(0.1*ParamPhysicalModel.P0_dBm)*1e-3; 

        for numChannel = 1:ParamSig.numChannel  
            if ParamControl.resample_as_FFT_or_Not
                Tx_Time_Data_temp = resample(Tx_Time_Data_CWDM{numChannel},ParamSys.Analog_Rate_afterDeMux,ParamDAC.DAC_Rate);
                Tx_Freq_Data_temp = fftshift(fft(Tx_Time_Data_temp));
                length_temp = length(Tx_Time_Data_temp);
                zero_insertion = zeros(1,(ParamSys.Analog_Rate_beforeDeMux/ParamSys.Analog_Rate_afterDeMux-1)*length_temp/2);
                Tx_Freq_Data_temp = [zero_insertion,Tx_Freq_Data_temp,zero_insertion];
                Tx_Time_Data_temp = ifft(ifftshift(Tx_Freq_Data_temp));
            else
                Tx_Time_Data_temp = resample(Tx_Time_Data_CWDM{numChannel},ParamSys.Analog_Rate_beforeDeMux,ParamDAC.DAC_Rate);
            end
            Time_Vector = 1/ParamSys.Analog_Rate_beforeDeMux*(0:(length(Tx_Time_Data_temp)-1));
            Tx_Time_Data_CWDM_mux = Tx_Time_Data_CWDM_mux + sqrt(ParamPhysicalModel.P0)*pwr_normalization(Tx_Time_Data_temp).*exp(1j*2*pi*ParamLas.deltaF(numChannel).*Time_Vector);
        end

        %% split step Fourier Transform (SSFT)    
        % *********************************************************************
        ParamPhysicalModel.AccumulateDistance = 0;
        dZ = ParamFib.FiberLength/ParamFib.numsec_SSFT/2;
        for idx = 1:ParamFib.numsec_SSFT 
            if ParamControl.Add_NL_or_Not
                if ParamControl.Runge_Kutta_or_Not 
                    [Tx_Time_Data_CWDM_mux,ParamPhysicalModel.AccumulateDistance ] = Add_NL_singleNLS_Runge_Kutta(Tx_Time_Data_CWDM_mux,ParamFib,dZ,ParamPhysicalModel.AccumulateDistance ); 
                else
                    ParamFib.Leff = (1-exp(-ParamFib.Loss_alpha*dZ))/ParamFib.Loss_alpha;
                    NL_phase = -1i*ParamFib.gamma*ParamFib.Leff*...
                               abs(Tx_Time_Data_CWDM_mux).^2;
                    Tx_Time_Data_CWDM_mux = Tx_Time_Data_CWDM_mux.*exp(NL_phase).*exp(-ParamFib.Loss_alpha/2*dZ);
                end   
            end

            if ParamControl.Add_CD_or_Not
                if ParamControl.Add_CD_individually_or_Not == 1
                    for numChannel = 1:ParamSig.numChannel
                        Time_Vector = 1/ParamSys.Analog_Rate_beforeDeMux*(0:(length(Tx_Time_Data_CWDM_mux)-1));
                        Tx_Time_Data_temp = Tx_Time_Data_CWDM_mux.*exp(-1j*2*pi*ParamLas.deltaF(numChannel).*Time_Vector);
                        if ParamControl.resample_as_FFT_or_Not
                            Tx_Freq_Data_temp = fftshift(fft(Tx_Time_Data_temp));
                            length_temp = length(Tx_Time_Data_temp);
                            num_zero_discard = (ParamSys.Analog_Rate_beforeDeMux/ParamSys.Analog_Rate_afterDeMux-1)/...
                                (ParamSys.Analog_Rate_beforeDeMux/ParamSys.Analog_Rate_afterDeMux)*length_temp/2;               
                            Tx_Freq_Data_temp(1:num_zero_discard) = [];
                            Tx_Freq_Data_temp(end-num_zero_discard+1:end) = [];
                            Tx_Time_Data_temp = ifft(ifftshift(Tx_Freq_Data_temp));
                            Tx_Time_Data_CWDM{numChannel} = Tx_Time_Data_temp;
                        else
                            Tx_Time_Data_CWDM{numChannel} = resample(Tx_Time_Data_temp,ParamSys.Analog_Rate_afterDeMux,ParamSys.Analog_Rate_beforeDeMux);    
                        end
                    end

                    DTime = 1/ParamSys.Analog_Rate_afterDeMux;
                    for numChannel = 1:ParamSig.numChannel
                        % consider vg mismatch
                        Tx_Time_Data_CWDM{numChannel} = Add_CD_TOD(Tx_Time_Data_CWDM{numChannel},...
                            DTime, dZ*2, ParamFib.Beta1(numChannel),...
                            ParamFib.Beta2(numChannel),ParamFib.Beta3(numChannel));
                        % consider phase mismatch
                        Tx_Time_Data_CWDM{numChannel} = Tx_Time_Data_CWDM{numChannel}.*exp(-1i*ParamFib.Beta0(numChannel)*dZ*2);
                    end

                    Tx_Time_Data_CWDM_mux = 0;
                    for numChannel = 1:ParamSig.numChannel
                        if ParamControl.resample_as_FFT_or_Not
                            Tx_Time_Data_temp = Tx_Time_Data_CWDM{numChannel};
                            Tx_Freq_Data_temp = fftshift(fft(Tx_Time_Data_temp));
                            length_temp = length(Tx_Time_Data_temp);
                            zero_insertion = zeros(1,(ParamSys.Analog_Rate_beforeDeMux/ParamSys.Analog_Rate_afterDeMux-1)*length_temp/2);
                            Tx_Freq_Data_temp = [zero_insertion,Tx_Freq_Data_temp,zero_insertion];
                            Tx_Time_Data_temp = ifft(ifftshift(Tx_Freq_Data_temp));
                        else
                            Tx_Time_Data_temp = resample(Tx_Time_Data_CWDM{numChannel},ParamSys.Analog_Rate_beforeDeMux,ParamSys.Analog_Rate_afterDeMux);
                        end
                        Time_Vector = 1/ParamSys.Analog_Rate_beforeDeMux*(0:(length(Tx_Time_Data_temp)-1));
                        Tx_Time_Data_CWDM_mux = Tx_Time_Data_CWDM_mux + Tx_Time_Data_temp.*exp(1j*2*pi*ParamLas.deltaF(numChannel).*Time_Vector);
                    end
                else
                    DTime = 1/ParamSys.Analog_Rate_beforeDeMux;
                    Tx_Time_Data_CWDM_mux = Add_CD_5th_order(Tx_Time_Data_CWDM_mux,...
                            DTime, dZ*2, ParamFib.Beta1_ref,...
                            ParamFib.Beta2_ref,ParamFib.Beta3_ref,ParamFib.Beta4_ref,ParamFib.Beta5_ref);

                end
            end

            if ParamControl.Add_NL_or_Not
                if ParamControl.Runge_Kutta_or_Not 
                    [Tx_Time_Data_CWDM_mux,ParamPhysicalModel.AccumulateDistance ] = Add_NL_singleNLS_Runge_Kutta(Tx_Time_Data_CWDM_mux,ParamFib,dZ,ParamPhysicalModel.AccumulateDistance ); 
                else
                    ParamFib.Leff = (1-exp(-ParamFib.Loss_alpha*dZ))/ParamFib.Loss_alpha;
                    NL_phase = -1i*ParamFib.gamma*ParamFib.Leff*...
                               abs(Tx_Time_Data_CWDM_mux).^2;
                    Tx_Time_Data_CWDM_mux = Tx_Time_Data_CWDM_mux.*exp(NL_phase).*exp(-ParamFib.Loss_alpha/2*dZ);
                end            
            end
        end

        if ParamControl.Plot_Muxed_Spectrum_or_Not 
            figure;
            Freq_Vector = linspace(-ParamSys.Analog_Rate_beforeDeMux/2,ParamSys.Analog_Rate_beforeDeMux/2,length(Tx_Time_Data_CWDM_mux));
            plot(Freq_Vector/1e9,10*log10(abs(FT(pwr_normalization(Tx_Time_Data_CWDM_mux),1/ParamSys.Analog_Rate_beforeDeMux))),'LineWidth',2);
            title('CWDM spectrum');
        end
        if ParamControl.save_Waveform_or_Not
            save('Temp_Param.mat','Tx_Sequence_CWDM','ParamPhysicalModel','ParamPS');
            save('Temp_Waveform.mat','Tx_Time_Data_CWDM_mux','-v7.3');
        end
    else
        ParamPhysicalModel.P0_dBm = ParamPhysicalModel.LchPwr_dBm - ParamChan.excess_loss_dB/2; 
        ParamPhysicalModel.P0 = 10^(0.1*ParamPhysicalModel.P0_dBm)*1e-3;
        %% Generate Analog waveform
        for numChannel = 1:ParamSig.numChannel
            Tx_Time_Data_temp = resample(Tx_Time_Data_CWDM{numChannel},ParamSys.Analog_Rate_afterDeMux,ParamDAC.DAC_Rate);
            Tx_Time_Data_CWDM{numChannel} = pwr_normalization(Tx_Time_Data_temp)*sqrt(ParamPhysicalModel.P0);
        end
        %% SSFT
        DTime = 1/ParamSys.Analog_Rate_afterDeMux;
        ParamPhysicalModel.AccumulateDistance = 0;
        dZ = ParamFib.FiberLength/ParamFib.numsec_SSFT/2;
        for idx = 1:ParamFib.numsec_SSFT 
            if ParamControl.Add_NL_or_Not  
                switch ParamControl.RK4_version
                    case 1
                        [Tx_Time_Data_CWDM,ParamPhysicalModel.AccumulateDistance ] = Add_NL_coupledNLS_Runge_Kutta(Tx_Time_Data_CWDM,ParamFib,dZ,ParamPhysicalModel.AccumulateDistance ); 
                    case 2
                        [Tx_Time_Data_CWDM,ParamPhysicalModel.AccumulateDistance ] = Add_NL_coupledNLS_Runge_Kutta_v2(Tx_Time_Data_CWDM,ParamFib,dZ,ParamPhysicalModel.AccumulateDistance ); 
                end
            end
            if ParamControl.Add_CD_or_Not          
                for numChannel = 1:ParamSig.numChannel
                    Tx_Time_Data_CWDM{numChannel} = Add_CD_TOD(Tx_Time_Data_CWDM{numChannel},...
                        DTime, dZ*2, ParamFib.Beta1(numChannel),...
                        ParamFib.Beta2(numChannel),ParamFib.Beta3(numChannel));
                    if ParamControl.RK4_version == 2
                        Tx_Time_Data_CWDM{numChannel} = Tx_Time_Data_CWDM{numChannel}.*exp(-1i*ParamFib.Beta0(numChannel)*dZ*2);
                    end
                end 
            end
            if ParamControl.Add_NL_or_Not
                switch ParamControl.RK4_version
                    case 1
                        [Tx_Time_Data_CWDM,ParamPhysicalModel.AccumulateDistance ] = Add_NL_coupledNLS_Runge_Kutta(Tx_Time_Data_CWDM,ParamFib,dZ,ParamPhysicalModel.AccumulateDistance ); 
                    case 2
                        [Tx_Time_Data_CWDM,ParamPhysicalModel.AccumulateDistance ] = Add_NL_coupledNLS_Runge_Kutta_v2(Tx_Time_Data_CWDM,ParamFib,dZ,ParamPhysicalModel.AccumulateDistance ); 
                end
            end
        end
        if ParamControl.save_Waveform_or_Not
            save('Temp_Param.mat','Tx_Sequence_CWDM','ParamPhysicalModel','ParamPS');
            save('Temp_Waveform.mat','Tx_Time_Data_CWDM','-v7.3');
        end
    end


    % *************************************************************************
    % *************************************************************************
    % *************************************************************************
else
    load('Temp_Waveform.mat');
    load('Temp_Param.mat');
end
%% Rx side
ParamPhysicalModel.Received_Pwr_dB = ParamPhysicalModel.LchPwr_dBm - ParamChan.Fiber_attenuation - ParamChan.excess_loss_dB; 
ParamPhysicalModel.BER_CWDM = [];


for numChannel = 1:ParamSig.numChannel
    if ParamControl.Coupled_NLS_or_Not == 0
        %% DeMux
        Time_Vector = 1/ParamSys.Analog_Rate_beforeDeMux*(0:(length(Tx_Time_Data_CWDM_mux)-1));
        Tx_Time_Data = Tx_Time_Data_CWDM_mux.*exp(-1j*2*pi*ParamLas.deltaF(numChannel).*Time_Vector);
        Tx_Time_Data = resample(Tx_Time_Data,ParamSys.Analog_Rate_afterDeMux,ParamSys.Analog_Rate_beforeDeMux);
    else
        Tx_Time_Data = Tx_Time_Data_CWDM{numChannel};
    end
    %% PD model
    % *********************************************************************
    Tx_Time_Data = pwr_normalization(Tx_Time_Data);
    Rx_I = Tx_Time_Data.*conj(Tx_Time_Data);
    % PD noise
    if ParamControl.Add_PD_noise
        Received_Pwr = 10^(0.1*ParamPhysicalModel.Received_Pwr_dB)*1e-3;%[W]
        switch ParamControl.PD_case    
            case 1 % PIN
                delta_f = ParamPD.BW_pin;
                SNR_PD_Linear = (ParamPD.R_pin*Received_Pwr)^2/...
                                (2*ParamGen.q_electron*(ParamPD.R_pin*Received_Pwr+ParamPD.Id_pin)*delta_f+...
                                 4*ParamGen.kB*ParamGen.Temperature/ParamPD.RL_pin*ParamPD.Fn_pin*delta_f);
                Rx_I = Rx_I*Received_Pwr*ParamPD.R_pin; 
                
                ParamPhysicalModel.SNR_PD_dB = 10*log10(SNR_PD_Linear);
                Rx_I = resample(Rx_I,delta_f*2,ParamSys.Analog_Rate_afterDeMux);
                IFFT_bin_length = 2*delta_f/ParamSig.Baud_Rate;
                Rx_I = Add_AWGN_real(Rx_I,ParamPhysicalModel.SNR_PD_dB,IFFT_bin_length,1);
                Rx_I = resample(Rx_I,ParamSys.Analog_Rate_afterDeMux,delta_f*2);
                
                Rx_V = Rx_I*ParamPD.RL_pin;
                
            case 2 % APD
                delta_f = ParamPD.BW_apd;
                ParamPD.FA_apd = ParamPD.kA_apd*ParamPD.M_apd+(1-ParamPD.kA_apd)*(2-1/ParamPD.M_apd);
                SNR_PD_Linear = (ParamPD.M_apd*ParamPD.R_apd*Received_Pwr)^2/...
                                (2*ParamGen.q_electron*ParamPD.M_apd^2*ParamPD.FA_apd*(ParamPD.R_apd*Received_Pwr+ParamPD.Id_apd)*delta_f+...
                                 4*ParamGen.kB*ParamGen.Temperature/ParamPD.RL_apd*ParamPD.Fn_apd*delta_f);
                Rx_I = Rx_I*Received_Pwr*ParamPD.R_apd*ParamPD.M_apd; 
                
                ParamPhysicalModel.SNR_PD_dB = 10*log10(SNR_PD_Linear);
                Rx_I = resample(Rx_I,delta_f*2,ParamSys.Analog_Rate_afterDeMux);
                IFFT_bin_length = 2*delta_f/ParamSig.Baud_Rate;
                Rx_I = Add_AWGN_real(Rx_I,ParamPhysicalModel.SNR_PD_dB,IFFT_bin_length,1);
                Rx_I = resample(Rx_I,ParamSys.Analog_Rate_afterDeMux,delta_f*2);
                
                Rx_V = Rx_I*ParamPD.RL_apd;
                
                
            case 3 % PIN-TIA NEP model
                Rx_V = Rx_I*Received_Pwr*ParamPD.CG; 
                PD_noise_Vrms = ParamPD.CG*ParamPD.NEP; 
                ParamPhysicalModel.SNR_PD_dB = 20*log10(rms(Rx_V)/PD_noise_Vrms);
                RL = 50;
                PD_noise = wgn(1,length(Rx_V),10*log10(PD_noise_Vrms^2/RL),RL,'dBW','real');
                Rx_V = Rx_V + PD_noise;               
        end
        ParamPD.SNR_PD_dB = ParamPhysicalModel.SNR_PD_dB;
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
        Rx_V = lpfilt(Rx_V,ParamSys.Analog_Rate_afterDeMux,PD_BW,'gaussian',PD_BW_ord);
    end
    
    %% ADC model
    % *********************************************************************
    % ADC noise (brick wall filter is also considered)
    if ParamControl.Add_ADC_thermal_noise
        Rx_V = resample(Rx_V,ParamADC.ADC_BW*2,ParamSys.Analog_Rate_afterDeMux);
        ADC_noise_power_LN = 10^(ParamADC.ADC_noise_power_dBm/10)/12.5e9*ParamADC.ADC_BW;
        noise = wgn(1,length(Rx_V),10*log10(ADC_noise_power_LN),50,'dBm','real');
        %ADC_noise_rms = rms(noise-mean(noise));
        %ADC_noise_Vpp = max(noise)-min(noise);
        Rx_V = Rx_V + noise;
        ParamADC.SNR_ADC = 10*log10(get_pwr(Rx_V)/get_pwr(noise));
        Rx_V = resample(Rx_V,ParamSys.Analog_Rate_afterDeMux,ParamADC.ADC_BW*2);   
    end
    
    % ADC detection
    Rx_V = resample(Rx_V,ParamADC.ADC_Rate,ParamSys.Analog_Rate_afterDeMux);
    
    ParamPhysicalModel.RV_rms_at_ADC = rms(Rx_V-mean(Rx_V));
    ParamPhysicalModel.RVpp_at_ADC = max(Rx_V)-min(Rx_V);
    
    % Rx_V = circshift(Rx_V,[0,2e3+2]);
    Rx_V = Quantization(real(Rx_V),ParamADC.qnbit_ADC);
    % Rx Side
    
    [ParamPhysicalModel.BER_SCMavg, ParamPhysicalModel.BER_list, ParamPhysicalModel.SNR_list] = zp_SP_SCM_DSB_Rx(ParamControl,ParamRxDSP,ParamSig,ParamADC,ParamPS,ParamFib,ParamMod,ParamPD,Rx_V,Tx_Sequence_CWDM{numChannel});  
    
    ParamPhysicalModel.BER_CWDM = [ParamPhysicalModel.BER_CWDM,ParamPhysicalModel.BER_SCMavg];
end

    ParamPhysicalModel.BER_CWDM_avg = mean(ParamPhysicalModel.BER_CWDM);

toc
% ParamGen
% ParamSig
% ParamDAC
% ParamLas
% ParamMod
% ParamFib
% ParamChan
% ParamVSB
% ParamPD
% ParamADC
% ParamSys
% ParamRxDSP
%% General
ParamGen.NoiseReferenceBW = 12.5e9; % 0.1 nm
ParamGen.c0 = 3e8;
ParamGen.Temperature = 273 + 23; %[K]
ParamGen.kB = 1.38064852e-23; %[m^2*kg*s^-2*K^-1] Boltzmann constant
ParamGen.q_electron = 1.60217662e-19;%[coulomb]

%% DAC init
ParamDAC.maxSamplesToLoad = 2^19;
ParamDAC.freq_resolution = 0.1e9;
switch ParamControl.FEC_option
    case 1
        ParamDAC.DAC_Rate = 64e9;
    case 2
        ParamDAC.DAC_Rate = 75e9;
end
%% Signal
switch ParamControl.PAM_or_QAM
    case 1
        if ParamControl.SCM_or_Not==0 || ParamControl.Curve_Fitting_or_Not
            switch ParamControl.FEC_option
                case 1
                    ParamSig.SC_Baud_Rate = 56e9;
                case 2
                    ParamSig.SC_Baud_Rate = 60e9;
            end
            ParamSig.SC = 2;
            ParamSig.SC_power = 1;

        else
            % currently PAM SCM still under development
        end
        ParamSig.GuardBand = 0; % just to avoid bugs in the code
        ParamSig.Baud_Rate = sum(ParamSig.SC_Baud_Rate);

        ParamSig.roll_off = 0.5;
    case 2
        if ParamControl.SCM_or_Not==0 || ParamControl.Curve_Fitting_or_Not
            switch ParamControl.FEC_option
                case 1
                    ParamSig.SC_Baud_Rate = 28e9;
                case 2
                    ParamSig.SC_Baud_Rate = 30e9;
            end
            ParamSig.SC = 4;
            ParamSig.SC_power = 1;

        else
            ParamSig.SC_Baud_Rate = [8,20]*1e9;
            ParamSig.SC = [3,4];
            ParamSig.SC_power = ([1.2,2.5]);
        end
        ParamSig.roll_off = 0.1;
        ParamSig.Baud_Rate = sum(ParamSig.SC_Baud_Rate);
        switch ParamControl.FEC_option
             
            case 1
                % range 247- 265  0.35 -> 252
                ParamSig.Ncircshift = 256; % 256
            
                ParamSig.GuardBand = ParamDAC.DAC_Rate/1024*...
                    ParamSig.Ncircshift-...
                    ParamSig.Baud_Rate*(1+ParamSig.roll_off)/2;
                
            case 2
                ParamSig.GuardBand = 0.140625e9;
        end
        if ParamControl.Curve_Fitting_or_Not
            ParamSig.GuardBand = 3e9;
        end
        
end
%% DAC (cont'd)


if ParamControl.Curve_Fitting_or_Not 
    ParamDAC.DAC_Rate = 88e9;
end
ParamDAC.qnbit_DAC = 8;
ParamDAC.clipping_Prob = 1e-3; % 1e-2

switch ParamControl.DAC_LPF_option
    case 1
        ParamDAC.DAC_LPF_BW = 12e9;
        ParamDAC.DAC_LPF_order = 0.6;
        ParamDAC.DAC_LPF_type = 'gaussian';
    case 2
        load('equalizer_D0IN_SHF39358_totkuVcable_DCA_electrical_head_88GSps_v2.mat')
%         load('G:\Project\active\IMDD_VSB_CDC\ZPOOLA\equalziers\equalizer_D0IN_SHF39358_totkuVcable_DCA_electrical_head_88GSps_v2.mat')

        eq_temp=TxRxWaveformat88Gsps.besth.';
        eq_temp_freq = fftshift(fft(eq_temp));
        BW_Hz = ParamSig.Baud_Rate*(1+ParamSig.roll_off)+ParamSig.GuardBand;
        N = length(eq_temp);
        if mod(N,2) == 0
            f = 88e9/N * (( -N/2 : N/2-1 ));
        else
            f = 88e9/N * ( -(N-1)/2 : (N-1)/2 );
        end
        indexPosChoose = f > BW_Hz;
        indexNegChoose = f < -BW_Hz;
        eq_temp_freq(indexPosChoose) = 0;
        eq_temp_freq(indexNegChoose) = 0;
        dac_lpf_temp_freq = zeros(size(eq_temp_freq));
        dac_lpf_temp_freq(f<=BW_Hz & f>=-BW_Hz) = 1./eq_temp_freq(f<=BW_Hz & f>=-BW_Hz);
        
        figure; plot(20*log10(abs(dac_lpf_temp_freq/max(dac_lpf_temp_freq))));
        ParamDAC.DAC_response_FIR = real(ifft(ifftshift(dac_lpf_temp_freq)));
        ParamDAC.Preemphasis_FIR = real(ifft(ifftshift(eq_temp_freq)));
        
end

ParamDAC.DAC_SNR= 30; % DAC noise
ParamDAC.DAC_SNR = ParamDAC.DAC_SNR + ...
                   10*log10(88e9/ParamDAC.DAC_Rate);
% **************
% % ParamADC.DAC_noise_power_dBm = -48.5; % per 0.1 nm 
% %                                       % consider Rx_V, 50 Ohm impedance

                                      
%% laser
switch ParamControl.Laser_case
    case 1
        ParamLas.laser_power_dBm = 10;
        ParamLas.Laser_Linewidth = 5e6;
    case 2
        ParamLas.laser_power_dBm = 13;
        ParamLas.Laser_Linewidth = 2e6;
    case 3
        if ParamControl.VSB_or_Not
            ParamLas.laser_power_dBm = 17.5;
        else
            ParamLas.laser_power_dBm = 14.5;
        end
        ParamLas.Laser_Linewidth = 1e6;
end
if ParamControl.Curve_Fitting_or_Not
    ParamLas.Laser_Linewidth = 200e3;
end

ParamLas.slope = 0.137; %10^(1.6)/(300-15);
ParamLas.current = 10^(ParamLas.laser_power_dBm/10)/ParamLas.slope+40;
ParamLas.Voltage = 2.7;
ParamLas.Laser_OSNR = 39; % laser OSNR 10*log10(10^(0.1*(-140))*50e9/4)= -39
if ParamControl.RxPreAmp_or_Not == 1
    ParamLas.Laser_OSNR =  ParamLas.Laser_OSNR - 5;
end

ParamLas.Lambda_ref = 1.55e-6; % C band
ParamLas.Lambda0 = 1.55e-6; % C band
ParamLas.omega0 = 2*pi*ParamGen.c0/ParamLas.Lambda0;
ParamLas.omega_ref = 2*pi*ParamGen.c0/ParamLas.Lambda_ref;

%% modulator
ParamMod.Modulator_Loss_dB = 5; % single modulator
ParamMod.Y_branch_Loss_dB = 0.5;
ParamMod.MMI_Loss_dB = 0.7;
ParamMod.Modulator_LPF_BW = 44e9;
ParamMod.Modulator_LPF_order = 6;
ParamMod.SSC_Loss_dB = 2;
ParamMod.carrier_path_loss_dB = 0.1;
%% Fiber
ParamFib.FiberLength = 40e3;
if ParamControl.Curve_Fitting_or_Not
    ParamFib.FiberLength = 80e3;
end

ParamFib.Loss_alpha_dB = 0.2; % [dB/km]
ParamFib.Loss_alpha = -log(10^(-0.1*ParamFib.Loss_alpha_dB/1e3));

zp_Fib_dispersion_Param_Init;

ParamFib.numsec_SSFT = 4;
%% Channel
ParamChan.SNR_penalty = 0;%[dB]

ParamChan.SNR = ParamLas.Laser_OSNR  + ...
                10*log10(ParamGen.NoiseReferenceBW/ParamSig.Baud_Rate)-...
                ParamChan.SNR_penalty;
            
ParamChan.coupler_loss_dB = 0.3;
ParamChan.extra_loss_dB = 0;
ParamChan.Fiber_attenuation = ParamFib.FiberLength*1e-3*ParamFib.Loss_alpha_dB;
            
%% VSB 
if ParamControl.VSB_or_Not
    switch ParamControl.OBPF_option 
        case 1

            ParamVSB.Opt_Flt_Suppression_dB = 50;
            ParamVSB.Opt_Flt_Slope_dBper10GHz = 1000; 
            ParamVSB.Opt_Flt_offset = 0e9; 
            % ParamVSB.Opt_Flt_offset =ParamSig.GuardBand;
            ParamVSB.Opt_Flt_drift = 0e9; 
            ParamVSB.Opt_Flt_ILoss_dB = 2.5 + ParamChan.coupler_loss_dB*2;
        case 5
            ParamVSB.Opt_Flt_offset = 3e9; 
            % ParamVSB.Opt_Flt_offset =ParamSig.GuardBand;
            ParamVSB.Opt_Flt_drift = 0e9; 
            ParamVSB.Opt_Flt_ILoss_dB = 0.5 + ParamChan.coupler_loss_dB*2;
    end
else
    ParamVSB = 0;
end
%% PD

switch ParamControl.PD_case
    case 1
        % 1 PIN    
        ParamPD.R_pin = 0.65; %[A/W]
        ParamPD.Id_pin = 5e-9; % [A]
        ParamPD.RL_pin = 50;%[Ohm] 
        ParamPD.Fn_pin = 1; % noise figure of electrical amplifier 
        ParamPD.BW_pin = 50e9; %[Hz]
        ParamPD.BW_PIN_order = 0.5;
    case 2
        % 2 APD
        ParamPD.R_apd = 0.6;
        ParamPD.Id_apd = 270e-9;
        ParamPD.RL_apd = 50;
        ParamPD.Fn_apd = 1;
        ParamPD.BW_apd = 36e9;
        ParamPD.M_apd = 10;
        ParamPD.kA_apd = 0.1;
        ParamPD.BW_APD_order = 1.08;
    case 3
        % 3 PIN-TIA
        ParamPD.CG = 724; % V/W
        ParamPD.BW_PIN_TIA = 37e9;
        ParamPD.BW_PIN_TIA_order = 2.3;
        ParamPD.NEP = 3.4e-6; % W
end

%% ADC(RTO)

switch ParamControl.FEC_option
    case 1
        ParamADC.ADC_Rate = 64e9;
    case 2 
        ParamADC.ADC_Rate = 75e9;
end

if ParamControl.Curve_Fitting_or_Not 
    ParamADC.ADC_Rate = 160e9;
end

ParamADC.qnbit_ADC = 8;


ParamADC.ADC_BW = ParamADC.ADC_Rate/2;


ParamADC.ADC_noise_power_dBm = -48.5-0; % per 0.1 nm  % -149.47 dBm/Hz
                                      % consider Rx_V, 50 Ohm impedance

          
%% other system parameter
ParamSys.Analog_Rate = 200e9;

if ParamControl.Curve_Fitting_or_Not || ParamControl.VSB_or_Not 
    ParamSys.Path_Mismatch = 100e-12;%[s]
else
    ParamSys.Path_Mismatch = 0e-12;%[s]
end

if  ParamControl.Curve_Fitting_or_Not ||  ParamControl.CSPR_tuning_case == 1 ...
        || ParamControl.RxPreAmp_or_Not    
    ParamSys.Target_RxPwr_dBm = 7;
end

% CSPR, Vbias, Vpp, Vpi
switch ParamControl.CSPR_tuning_case
    case 1
        ParamSys.CSPR_dB = 12; % for linear model CSPR tuning
    case 2
        ParamSys.Carrier_path_pwr_ratio = 0.17;
        if ParamControl.VSB_or_Not
            ParamSys.Carrier_path_pwr_ratio = 0.34;
        end
        ParamSys.Vbias_over_Vpi = 0;
    case 3
        if ParamControl.PAM_or_QAM == 1
            ParamSys.Vbias_over_Vpi = 0.22;
        else     
            ParamSys.Vbias_over_Vpi = 0.44;
        end
end

ParamSys.Vpi = 1.3;
ParamSys.Vpp_over_Vpi = 1.4;

% if ParamControl.VSB_or_Not == 0
%     ParamSys.Vbias_over_Vpi = ParamSys.Vbias_over_Vpi*sqrt(2);
% end

% Not that when using KK we also nee8
% ParamSys.Vbias_over_Vpi > ParamSys.Vpp_over_Vpi/2;



if ParamControl.CSPR_tuning_case  == 3
    if ParamSys.Vpp_over_Vpi/2 > 1-ParamSys.Vbias_over_Vpi
        ParamSys.Vpp_over_Vpi = 2*(1-ParamSys.Vbias_over_Vpi);
        disp('V drive should not go to the second slope !');
    end
else
    ParamSys.Vbias_over_Vpi = 0;
end


%% RxDSP
if ParamControl.Digital_Resample_Before_KK_or_Not 
    switch ParamControl.FEC_option
        case 1
            if ParamControl.VSB_or_Not
                ParamRxDSP.KKoverSamp = 80/28;
            else
                ParamRxDSP.KKoverSamp = 64/28;
            end
        case 2
            if ParamControl.KK_option == 0
                ParamRxDSP.KKoverSamp = 90/30;
            else
                ParamRxDSP.KKoverSamp = 75/30;
            end
    end
else
    ParamRxDSP.KKoverSamp = (ParamADC.ADC_Rate/1e9)/(ParamSig.Baud_Rate/1e9);
end
ParamRxDSP.nb = 8;
ParamRxDSP.f_clk = 500e6; % 500 MHz
ParamRxDSP.hilbert_tap = 0; % =0 no overlap at all
ParamRxDSP.numKKiter = 1; % for KK option 2,3 and 4 
ParamRxDSP.FAiterKK_tap = 0;
if ParamControl.VSB_or_Not
    ParamRxDSP.CD_tap = 35; 
else
    ParamRxDSP.CD_tap = 28; 
end
ParamRxDSP.Nfft = 1024;
ParamRxDSP.FFT_size_ratio = 1; 
if ParamControl.VSB_or_Not
    ParamRxDSP.LMS_linear_tap = 3; 
else
    ParamRxDSP.LMS_linear_tap = 3; 
end

ParamRxDSP.LMS_DD_step = 5e-4;
ParamRxDSP.LMS_Train_step = 1e-3;
ParamRxDSP.Train_Sequence_Length = 2^14;
ParamRxDSP.Head = 1000;

if ParamControl.Curve_Fitting_or_Not
    ParamRxDSP.LMS_linear_tap = 51;
end


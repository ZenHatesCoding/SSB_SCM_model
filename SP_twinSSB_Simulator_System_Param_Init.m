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

%% Signal
if ~ParamControl.SCM_or_Not || ParamControl.Curve_Fitting_or_Not
    switch ParamControl.FEC_option
        case 1
            ParamSig.SC_Baud_Rate = 28e9;
        case 2
            ParamSig.SC_Baud_Rate = 30e9;
    end
    ParamSig.SC = 4;
    ParamSig.SC_power = 1;
   
else
    ParamSig.SC_Baud_Rate = [8,8,8,8]*1e9;
    ParamSig.SC = [4,4,4,4];
    ParamSig.SC_power = ([1,1,1,1]);
end

ParamSig.Baud_Rate = sum(ParamSig.SC_Baud_Rate);
switch ParamControl.FEC_option
    case 1
        ParamSig.GuardBand = 0.35e9;
    case 2
        ParamSig.GuardBand = 0.140625e9;
end
if ParamControl.Curve_Fitting_or_Not
    ParamSig.GuardBand = 3e9;
end
ParamSig.roll_off = 0.1;

%% DAC
ParamDAC.maxSamplesToLoad = 2^18;
ParamDAC.freq_resolution = 0.1e9;
switch ParamControl.FEC_option
    case 1
        ParamDAC.DAC_Rate = 64e9;
    case 2
        ParamDAC.DAC_Rate = 75e9;
end

if ParamControl.Curve_Fitting_or_Not 
    ParamDAC.DAC_Rate = 88e9;
end
ParamDAC.qnbit_DAC = 8;
ParamDAC.clipping_Prob = 10e-3;

switch ParamControl.DAC_LPF_option
    case 1
        ParamDAC.DAC_LPF_BW = 23e9;
        ParamDAC.DAC_LPF_order = 1;
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
        ParamDAC.DAC_response_FIR = real(ifft(ifftshift(dac_lpf_temp_freq)));
        ParamDAC.Preemphasis_FIR = real(ifft(ifftshift(eq_temp_freq)));
        
end

ParamDAC.DAC_SNR= 32; % DAC noise
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
        ParamLas.laser_power_dBm = 16;
        ParamLas.Laser_Linewidth = 1e6;
end
if ParamControl.Curve_Fitting_or_Not
    ParamLas.Laser_Linewidth = 200e3;
end
ParamLas.Laser_OSNR = 40; % laser OSNR
ParamLas.Lambda_ref = 1.55e-6; % C band
ParamLas.Lambda0 = 1.55e-6; % C band
ParamLas.omega0 = 2*pi*ParamGen.c0/ParamLas.Lambda0;
ParamLas.omega_ref = 2*pi*ParamGen.c0/ParamLas.Lambda_ref;

%% modulator
ParamMod.Modulator_Loss_dB = 5; % single modulator

ParamMod.Modulator_Loss_dB = ParamMod.Modulator_Loss_dB - 10*log10(1/2);
% because the parent modulator always has pi/2 phase difference between
% the 2 arms

ParamMod.Modulator_LPF_BW = 30e9;
ParamMod.Modulator_LPF_order = 1.4;

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
                10*log10(ParamGen.NoiseReferenceBW/ParamSig.Baud_Rate/2)-...
                ParamChan.SNR_penalty;
            
ParamChan.excess_loss_dB = 2;
ParamChan.Fiber_attenuation = ParamFib.FiberLength*1e-3*ParamFib.Loss_alpha_dB;
            
%% Optical filter 
ParamVSB.Opt_Flt_Suppression_dB1 = 30;
ParamVSB.Opt_Flt_Slope_dBper10GHz1 = 190; % for case 1 linear filter 
ParamVSB.Opt_Flt_offset1 = 0e9;

ParamVSB.Opt_Flt_Suppression_dB2 = 30;
ParamVSB.Opt_Flt_Slope_dBper10GHz2 = 190; % for case 1 linear filter 
ParamVSB.Opt_Flt_offset2 = 0e9;
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
        if ParamControl.KK_option == 1
            ParamADC.ADC_Rate = 64e9;
        else
            ParamADC.ADC_Rate = 64e9;
        end

    case 2
        if ParamControl.KK_option == 1
            ParamADC.ADC_Rate = 75e9;
        else
            ParamADC.ADC_Rate = 75e9;
        end

end

if ParamControl.Curve_Fitting_or_Not 
    ParamADC.ADC_Rate = 160e9;
end

ParamADC.qnbit_ADC = 8;
ParamADC.ADC_BW = 63e9; % brick wall
if ParamADC.ADC_Rate/2 < ParamADC.ADC_BW
    ParamADC.ADC_BW = ParamADC.ADC_Rate/2;
end

ParamADC.ADC_noise_power_dBm = -48.5; % per 0.1 nm 
                                      % consider Rx_V, 50 Ohm impedance

          
%% other system parameter
ParamSys.Analog_Rate = 250e9;

ParamSys.Path_Mismatch = 0;

if  ParamControl.Curve_Fitting_or_Not ||  ParamControl.CSPR_tuning_case == 1 ...
        || ParamControl.RxPreAmp_or_Not    
    ParamSys.Target_RxPwr_dBm = -7;
end

% CSPR, Vbias, Vpp, Vpi
switch ParamControl.CSPR_tuning_case
    case 1
        ParamSys.CSPR_dB = 11; % for linear model CSPR tuning
        ParamSys.Vbias_over_Vpi = 0;
    case 2
        ParamSys.Carrier_path_pwr_ratio = 0.1;
        ParamSys.Vbias_over_Vpi = 0;
    case 3
        ParamSys.Vbias_over_Vpi = 0.2;
end

ParamSys.Vpp_over_Vpi = 0.33; 
% Not that when using KK we also nee8
% ParamSys.Vbias_over_Vpi > ParamSys.Vpp_over_Vpi/2;

ParamSys.Vpi = 6;
if ParamControl.CSPR_tuning_case  == 3
    if ParamSys.Vpp_over_Vpi/2 > 1-ParamSys.Vbias_over_Vpi
        ParamSys.Vpp_over_Vpi = 2*(1-ParamSys.Vbias_over_Vpi);
        disp('V drive should not go to the second slope !');
    end
else
    ParamSys.Vbias_over_Vpi = 0;
end


%% RxDSP

switch ParamControl.FEC_option
    case 1
        if ParamControl.KK_option == 0
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
ParamRxDSP.hilbert_tap = 1; % =1 no overlap at all
ParamRxDSP.numKKiter = 1; % for KK option 2 and 3 
ParamRxDSP.FFT_size_ratio = 1; 
ParamRxDSP.LMS_linear_tap = 25;
ParamRxDSP.LMS_DD_step = 1e-4;
ParamRxDSP.LMS_Train_step = 1e-3;
ParamRxDSP.Train_Sequence_Length = 2^14;
ParamRxDSP.Head = 1000;

if ParamControl.Curve_Fitting_or_Not
    ParamRxDSP.LMS_linear_tap = 51;
end


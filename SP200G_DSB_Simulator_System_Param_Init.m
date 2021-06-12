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
ParamGen.h = 6.626e-34; %[J*s]

%% Signal

switch ParamControl.FEC_option
    case 1
        ParamSig.SC_Baud_Rate = 56e9;
    case 2
        ParamSig.SC_Baud_Rate = 60e9;
    case 3
        ParamSig.SC_Baud_Rate = 92e9;
end
ParamSig.SC = 2.5;

ParamSig.Ninterleave = 2; % interestingly, 2 seems to be the best

ParamSig.SC_power = 1;   
ParamSig.Baud_Rate = sum(ParamSig.SC_Baud_Rate);
ParamSig.roll_off = 0.5;


ParamSig.numChannel = 1; % for DSB_PhysicalLayerSimulator_CWDM only
%% DAC
ParamDAC.maxSamplesToLoad = 2^18;
ParamDAC.freq_resolution = 0.1e9;
switch ParamControl.FEC_option
    case 1
        ParamDAC.DAC_Rate = 64e9;
    case 2
        ParamDAC.DAC_Rate = 75e9;
    case 3
        ParamDAC.DAC_Rate = 138e9;
end

ParamDAC.qnbit_DAC = 8;
ParamDAC.clipping_Prob = 1e-3;

switch ParamControl.DAC_LPF_option
    case 1
        ParamDAC.DAC_LPF_BW = 23e9;
        ParamDAC.DAC_LPF_order = 1;
        ParamDAC.DAC_LPF_type = 'gaussian';
    case 2
%         load('C:\Users\zxing3\Documents\GitHub\IMDD_VSB_CDC\ZPOOLA\equalziers\equalizer_D0IN_SHF39358_totkuVcable_DCA_electrical_head_88GSps_v2.mat')
        load('G:\Project\active\IMDD_VSB_CDC\ZPOOLA\equalziers\equalizer_D0IN_SHF39358_totkuVcable_DCA_electrical_head_88GSps_v2.mat')

        eq_temp=TxRxWaveformat88Gsps.besth.';
        eq_temp_freq = fftshift(fft(eq_temp));
        BW_Hz = ParamSig.Baud_Rate*(1+ParamSig.roll_off)/2;
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

                                      
%% laser
switch ParamControl.Laser_case
    case 1
        ParamLas.laser_power_dBm = 10;
        ParamLas.Laser_Linewidth = 5e6;
    case 2
        ParamLas.laser_power_dBm = 15;
        ParamLas.Laser_Linewidth = 2e6;
    case 3
        ParamLas.laser_power_dBm = 16;
        ParamLas.Laser_Linewidth = 1e6;
end

ParamLas.Laser_OSNR = 40; % laser OSNR
ParamLas.Lambda0 = 1.28e-6; % O band
ParamLas.Lambda_ref = ParamLas.Lambda0 ; % O band

ParamLas.omega0 = 2*pi*ParamGen.c0/ParamLas.Lambda0;
ParamLas.omega_ref = 2*pi*ParamGen.c0/ParamLas.Lambda_ref;

ParamLas.Photon_Energy = ParamLas.omega0/2/pi*ParamGen.h;
%% modulator
ParamMod.Modulator_Loss_dB = 5; % single modulator

ParamMod.Modulator_LPF_BW = 45e9;
ParamMod.Modulator_LPF_order = 1.4;

%% Fiber
ParamFib.FiberLength = 2e3;
ParamFib.Loss_alpha_dB = 0.35; % [dB/km]
ParamFib.Loss_alpha = -log(10^(-0.1*ParamFib.Loss_alpha_dB/1e3));

zp_Fib_dispersion_Param_Init;
ParamFib.numsec_SSFT = 4;

%% Channel
ParamChan.SNR_penalty = 0;%[dB]

ParamChan.SNR = ParamLas.Laser_OSNR  + ...
                10*log10(ParamGen.NoiseReferenceBW/ParamSig.Baud_Rate)-...
                ParamChan.SNR_penalty;
            
ParamChan.excess_loss_dB = 4;
ParamChan.Fiber_attenuation = ParamFib.FiberLength*1e-3*ParamFib.Loss_alpha_dB;
            

%% PD

switch ParamControl.PD_case
    case 1
        % 1 PIN    
        % ParamPD.R_pin = 0.65; %[A/W]
        ParamPD.R_pin = 0.55 * ParamGen.q_electron/ParamLas.Photon_Energy; %[A/W]
        ParamPD.Id_pin = 5e-9; % [A]
        ParamPD.RL_pin = 50;%[Ohm] 
        ParamPD.Fn_pin = 1; % noise figure of electrical amplifier 
        ParamPD.BW_pin = 50e9; %[Hz]
        ParamPD.BW_PIN_order = 0.5;
        
        ParamPD.Vbias = 2.8;
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
        ParamPD.R_pin = 0.55 * ParamGen.q_electron/ParamLas.Photon_Energy; %[A/W]

        ParamPD.CG = 724*ParamPD.R_pin/0.76; % V/W
        ParamPD.BW_PIN_TIA = 50e9;
        ParamPD.BW_PIN_TIA_order = 2.3;
        ParamPD.NEP = 3.4e-6*ParamPD.R_pin/0.71; % W % keep RMS noise to be the same
        ParamPD.Vbias = 5;
end

%% ADC(RTO)

switch ParamControl.FEC_option
    case 1
        ParamADC.ADC_Rate = 64e9;
    case 2
        ParamADC.ADC_Rate = 75e9;
    case 3
        ParamADC.ADC_Rate = 138e9;
end



ParamADC.qnbit_ADC = 8;
ParamADC.ADC_BW = 63e9; % brick wall
if ParamADC.ADC_Rate/2 < ParamADC.ADC_BW
    ParamADC.ADC_BW = ParamADC.ADC_Rate/2;
end

ParamADC.ADC_noise_power_dBm = -48.5; % per 0.1 nm 
                                      % consider Rx_V, 50 Ohm impedance

          
%% other system parameter
if ParamSig.numChannel == 1
    ParamSys.Analog_Rate = 250e9;
else
    ParamSys.Analog_Rate = 250e9*40;
end


% Vbias, Vpp, Vpi
ParamSys.Vbias_over_Vpi = 0.5;
ParamSys.Vpp_over_Vpi = 0.8; 
% Not that when using KK we also nee8
% ParamSys.Vbias_over_Vpi > ParamSys.Vpp_over_Vpi/2;
ParamSys.Vpi = 2.5;

%% RxDSP
ParamRxDSP.FFT_size_ratio = 1; 
ParamRxDSP.LMS_linear_tap = 5;
ParamRxDSP.LMS_DD_step =5e-4;
ParamRxDSP.LMS_Train_step = 1e-3;
ParamRxDSP.Train_Sequence_Length = 2^14;
ParamRxDSP.Head = 1024;



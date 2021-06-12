
%% signal
ParamControl.SCM_or_Not = 0;

%% Tx DSP
ParamControl.TxDSP_practical_implementation_or_Not =0;
ParamControl.Pulseshaping_or_Not = 1; % 1 using RRC/RC for pulse shaping 
                                      % 0 NRZ
    ParamControl.Use_RC_or_RRC = 0; % 0 RRC; 1 RC
ParamControl.CD_Pre_Compensation_or_Not = 0;
ParamControl.Preemphasis_or_Not = 1;
ParamControl.Clipping_or_Not = 0;

%% Physical layer
ParamControl.VSB_or_Not = 1; % related to how many dac channels to use
    ParamControl.OBPF_option = 2; % 0 rect, 1, linear, 2 practical, 3 SiPh
% DAC
ParamControl.DAC_LPF_or_Not = 1;
    ParamControl.DAC_LPF_option = 2; % 1 gaussian,2 practical
    % option 1 has some issue
ParamControl.Add_DAC_noise = 1;


% modulator
ParamControl.Modulator_Nonlinearity_or_Not = 0;
ParamControl.Modulator_LPF_or_Not = 1;

ParamControl.CSPR_tuning_case = 1; % 1 Linear model 2 Parallel Path 3 bias tuning

% laser
ParamControl.Add_PhaseNoise_or_Not = 1;
% Channel
ParamControl.Add_AWGN_or_Not = 1;
ParamControl.Measure_CSPR_or_Not = 1;
ParamControl.Add_CD_or_Not = 1;

% PD
ParamControl.Add_PD_noise = 1;
ParamControl.PD_LPF_or_Not = 1;
    ParamControl.PD_case = 1; % 1 PIN 2 APD 3 PIN_TIA 
    
% ADC
ParamControl.Add_ADC_thermal_noise = 1;

%% Rx DSP 
ParamControl.RxDSP_practical_implementation_or_Not = 0;
ParamControl.Use_DD_or_KK = 0; % DD = 1 KK = 0
    ParamControl.KK_option = 0; %0 normal KK, 1 KK w/o upsamp
ParamControl.CD_Compensation_or_Not = 1;
ParamControl.Synchronization_or_Not = 1;
ParamControl.Equalization_or_Not = 1;
    ParamControl.Update_Eqtap_or_Not = 1;
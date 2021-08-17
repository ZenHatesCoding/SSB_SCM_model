%% Block Control
ParamControl.PCA_Enable_or_Not = 1;

%% Plot
ParamControl.Plot_VSB_filter_or_Not = 0;
ParamControl.Plot_Spectrum_or_Not = 1;
ParamControl.Plot_CrossCorr_or_Not = 0;
ParamControl.Equalization_Plot_Conv_or_Not = 0;
ParamControl.Plot_Constellation_or_Not = 0;
ParamControl.Plot_Noise_or_Not = 0;
ParamControl.Plot_Error_or_Not = 0;
ParamControl.Plot_Power_Pie_or_Not = 0;

%% Experiment
ParamControl.new_load = 1;
ParamControl.new_capture  = 1;

ParamControl.Quantization_or_Not = 1;
%% signal
ParamControl.SCM_or_Not = 0;
ParamControl.PAM_or_QAM = 2; % 1 PAM 2 QAM
%% Tx DSP
ParamControl.TxDSP_practical_implementation_or_Not = 1;
if ParamControl.TxDSP_practical_implementation_or_Not == 1
    ParamControl.MergePreEmpwithUpconversion_or_Not = 1; % for QAM only
end
ParamControl.Pulseshaping_or_Not = 1; % 1 using RRC/RC for pulse shaping 
                                      % 0 NRZ
    ParamControl.Use_RC_or_RRC = 1; % 0 RRC; 1 RC
ParamControl.CD_Pre_Compensation_or_Not = 0;
ParamControl.Preemphasis_or_Not = 1;
ParamControl.Clipping_or_Not = 1;

%% Physical layer
ParamControl.VSB_or_Not = 0; % related to how many dac channels to use

    
% DAC
ParamControl.DAC_LPF_or_Not = 1;
    ParamControl.DAC_LPF_option = 1; % 1 gaussian,2 practical
    % option 1 has some issue
ParamControl.Add_DAC_noise = 1;


% modulator
ParamControl.Modulator_Nonlinearity_or_Not = 1;
if ParamControl.VSB_or_Not == 0
    ParamControl.MZM_option = 1; % 1 IQM 2 DDMZM
end

ParamControl.Modulator_LPF_or_Not = 1;

ParamControl.CSPR_tuning_case = 2; % 1 Linear model 2 Parallel Path 3 bias tuning

% laser
ParamControl.Laser_case = 3; % 1-10 dBm 2-13 dBm 3-16 dBm
ParamControl.Add_PhaseNoise_or_Not = 1;

% optical filter
if ParamControl.VSB_or_Not == 1
    ParamControl.OBPF_option = 5; % 0 rect, 1, linear, 2 practical, 3 SiPh
                                  % 5 lumentum
    ParamControl.LSB_or_RSB = 1; % 1 keep LSB 0 keep RSB
    ParamControl.VSB_at_Tx_or_Rx = 0; % 1 filter at Tx 0 at Rx
end
% Channel
ParamControl.Add_AWGN_beforeFiber_or_Not = 1;
ParamControl.Add_AWGN_afterFiber_or_Not = 0;
ParamControl.Measure_CSPR_or_Not = 1;
ParamControl.Add_CD_or_Not = 1;
    ParamControl.CD_order = 5;
ParamControl.Add_NL_or_Not = 1;

% PD
ParamControl.RxPreAmp_or_Not = 0;
ParamControl.Add_PD_noise = 1;
ParamControl.PD_LPF_or_Not = 1;
    ParamControl.PD_case = 3; % 1 PIN 2 APD 3 PIN_TIA 
    
% ADC
ParamControl.Add_ADC_thermal_noise = 0; % independent of input swing
ParamControl.Add_ADC_noise_distortion = 1; % defined as SNR
% For DSP
ParamControl.Add_internal_noise_or_Not = 1; % only for QAM now

%% Rx DSP 
ParamControl.Digital_Resample_Before_KK_or_Not = 1;
ParamControl.RxDSP_practical_implementation_or_Not = 1;

ParamControl.Use_DD_or_KK = 0; % DD = 1 KK = 0 (also includes SSBI mitigation)
    ParamControl.KK_option = 2; %0 normal KK, 
                                %1 KK w/o upsamp
                                %2 iterative KK
                                %3 iterKK for VSB 
                                %4 iterative SSBI mitigation

    if ParamControl.VSB_or_Not == 1
        ParamControl.KK_option = 3;
    end
    if ParamControl.KK_option == 0
        ParamControl.Hilbert_time_or_Freq = 0; % 1 time 0 freq 
    end
        ParamControl.CascadeKK_or_Not = 0; % stupid trick test

    ParamControl.KK_real_output_or_Not = 1;
                                         % always real output when using
                                         % iterKK
                                         
    ParamControl.FFT_sharing_or_Not = 1; % FFT sharing in resampling before KK
                                         % and hilbert transform
                                         % always share FFT when using
                                         % iterKK
    if ParamControl.KK_option == 2 || ParamControl.KK_option == 3
        ParamControl.FFT_sharing_or_Not = 1;
        ParamControl.KK_real_output_or_Not = 1;
    end
        ParamControl.Archive_or_SPPcom = 0; % for power consumption analysis only;
                                            % for HD-FEC only
                                            % for ufKK only
    if ParamControl.KK_real_output_or_Not ==1
        ParamControl.FFT_sharing_or_Not = 1;
        ParamControl.KK_real_output_subFFT_num = 1; % 1 or 2
    end
if ParamControl.RxDSP_practical_implementation_or_Not == 1
    ParamControl.MergeCDCwithdownconversion_or_Not = 1; % for QAM only
end

if ParamControl.RxDSP_practical_implementation_or_Not==0 ||...
        ParamControl.MergeCDCwithdownconversion_or_Not == 0
    ParamControl.CD_Compensation_or_Not = 1;
end


ParamControl.Synchronization_or_Not = 1;
ParamControl.Equalization_or_Not = 1;
    ParamControl.PLL_or_Not  = 0;
    ParamControl.Update_Eqtap_or_Not = 1;
    
%% CurveFitting
ParamControl.Curve_Fitting_or_Not = 0;
if ParamControl.Curve_Fitting_or_Not == 1
    ParamControl.CSPR_tuning_case = 1;
    ParamControl.TxDSP_practical_implementation_or_Not = 0;
    ParamControl.RxDSP_practical_implementation_or_Not = 0;
    ParamControl.Add_NL_or_Not = 0;
    ParamControl.RxPreAmp_or_Not = 1;
end
%% FEC option
ParamControl.FEC_option = 2;% 1 HD-FEC, 2 SD-FEC


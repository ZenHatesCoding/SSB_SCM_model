%% Block Control
ParamControl.PCA_Enable_or_Not = 0;

%% Plot
ParamControl.Plot_VSB_filter_or_Not = 0;
ParamControl.Plot_Spectrum_or_Not = 0;
ParamControl.Plot_CrossCorr_or_Not = 0;
ParamControl.Equalization_Plot_Conv_or_Not = 0;
ParamControl.Plot_Constellation_or_Not = 0;
ParamControl.Plot_Error_or_Not = 0;
ParamControl.Plot_Noise_or_Not = 0;
ParamControl.Plot_Power_Pie_or_Not = 0;

%% Experiment
ParamControl.new_load = 1;
ParamControl.new_capture  = 1;
ParamControl.Quantization_or_Not = 0;

%% signal
ParamControl.SCM_or_Not = 0;
ParamControl.PAM_or_QAM = 1; % 1 PAM 2 QAM

%% Tx DSP
ParamControl.TxDSP_practical_implementation_or_Not =0;
ParamControl.Pulseshaping_or_Not = 1; % 1 using RRC/RC for pulse shaping 
                                      % 0 NRZ
    ParamControl.Use_RC_or_RRC = 1; % 0 RRC; 1 RC
ParamControl.CD_Pre_Compensation_or_Not = 0;
ParamControl.Preemphasis_or_Not = 0;
ParamControl.Clipping_or_Not = 0;

%% Physical layer
% VSB
ParamControl.VSB_or_Not = 1; % related to how many dac channels to use
    ParamControl.OBPF_option = 1; % 0 rect, 1, linear, 2 practical, 3 SiPh
                                  % 4 ideal filter, stupid trick
ParamControl.VSB_at_Tx_or_Rx = 0; % 1 Tx 0 Rx
% DAC
ParamControl.DAC_LPF_or_Not = 0;
    ParamControl.DAC_LPF_option = 1; % 1 gaussian,2 practical
    % option 1 has some issue
ParamControl.Add_DAC_noise = 0;


% modulator
ParamControl.Modulator_Nonlinearity_or_Not = 0;
ParamControl.Modulator_LPF_or_Not = 0;

ParamControl.CSPR_tuning_case = 1; % 1 Linear model 2 Parallel Path 3 bias tuning

% laser
ParamControl.Laser_case = 2; % 1-10 dBm 2-13 dBm 3-16 dBm
ParamControl.Add_PhaseNoise_or_Not = 0;
% Channel
ParamControl.Add_AWGN_beforeFiber_or_Not = 0;
ParamControl.Add_AWGN_afterFiber_or_Not = 1;

ParamControl.Measure_CSPR_or_Not = 1;
ParamControl.Add_CD_or_Not = 1;
    ParamControl.CD_order = 2;
ParamControl.Add_NL_or_Not = 0;

% PD
ParamControl.RxPreAmp_or_Not = 1;
ParamControl.Add_PD_noise = 0;
ParamControl.PD_LPF_or_Not = 0;
    ParamControl.PD_case = 3; % 1 PIN 2 APD 3 PIN_TIA 
    
% ADC
ParamControl.Add_ADC_thermal_noise = 0;

%% Rx DSP 

ParamControl.RxDSP_practical_implementation_or_Not =0;
ParamControl.Use_DD_or_KK = 0; % DD = 1 KK = 0
    ParamControl.KK_option = 3; %0 normal KK, 1 KK w/o upsamp
                                %2 iterative KK 
                                %3 iterKK for VSB LN model (no practical impelmentation)
                                
        ParamControl.CascadeKK_or_Not = 0; % stupid trick test
    ParamControl.FFT_sharing_or_Not = 0; % FFT sharing in resampling before KK
                                         % and hilbert transform
                                         % always share FFT when using
                                         % iterKK
    ParamControl.KK_real_output_or_Not = 1;
                                         % always real output when using
                                         % iterKK
    if ParamControl.KK_option == 2
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
    
ParamControl.CD_Compensation_or_Not = 1;

    ParamControl.CDC_integratedwResamp_or_Not = 0; 
    % always set this to zero
    % this crap doesn't work,I keep this block just in case I wanna revist
    % this crap one day
ParamControl.Synchronization_or_Not = 1;
ParamControl.Equalization_or_Not = 1;
    ParamControl.PLL_or_Not  = 0;
    ParamControl.Update_Eqtap_or_Not = 1;
    
%% CurveFitting
ParamControl.Curve_Fitting_or_Not = 0;
%% FEC option
ParamControl.FEC_option = 1;% 1 HD-FEC, 2 SD-FEC


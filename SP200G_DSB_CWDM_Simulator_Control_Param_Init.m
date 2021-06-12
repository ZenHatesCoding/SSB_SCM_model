%% Block Control
ParamControl.PCA_Enable_or_Not = 0;
%% Plot Control
ParamControl.Plot_Spectrum_or_Not = 0;
    ParamControl.Plot_Muxed_Spectrum_or_Not = 0; % for single NLS only
ParamControl.Plot_CrossCorr_or_Not = 0;
ParamControl.Equalization_Plot_Conv_or_Not = 0;
ParamControl.Plot_Constellation_or_Not = 0;
ParamControl.Plot_Error_or_Not = 0;
ParamControl.Plot_Power_Pie_or_Not = 1;
%% Experiment

ParamControl.New_Simulation = 1;
    if ParamControl.New_Simulation == 1
        ParamControl.save_Waveform_or_Not = 1;
    end
%% signal
ParamControl.SCM_or_Not = 0;

%% Tx DSP
ParamControl.TxDSP_practical_implementation_or_Not =1;
ParamControl.Pulseshaping_or_Not = 1; % 1 using RRC/RC for pulse shaping 
                                      % 0 NRZ
    ParamControl.Use_RC_or_RRC = 1; % 0 RRC; 1 RC
ParamControl.CD_Pre_Compensation_or_Not = 0;
ParamControl.Preemphasis_or_Not = 1;
ParamControl.Clipping_or_Not = 1;

%%
ParamControl.IM_or_Not  = 1; % currently only consdier Intensity modulation
                             % no CD pre-compensation
%% Physical layer
% DAC
ParamControl.DAC_LPF_or_Not = 1;
    ParamControl.DAC_LPF_option = 1; % 1 gaussian,2 practical
    % option 1 has some issue
ParamControl.Add_DAC_noise = 1;


% modulator
ParamControl.Modulator_Nonlinearity_or_Not = 1;
ParamControl.Modulator_LPF_or_Not = 1;

% laser
ParamControl.Laser_case = 2; % 1-10 dBm 2-13 dBm 3-16 dBm
ParamControl.Add_PhaseNoise_or_Not = 1;
ParamControl.TEC_or_Not= 0;
% Channel
ParamControl.Add_AWGN_or_Not = 1;
ParamControl.Add_CD_or_Not = 1;

ParamControl.Add_NL_or_Not = 1;
ParamControl.Coupled_NLS_or_Not = 0;

if ParamControl.Coupled_NLS_or_Not == 0
    ParamControl.Add_CD_individually_or_Not = 0;

    ParamControl.resample_as_FFT_or_Not = 1;
    ParamControl.Runge_Kutta_or_Not = 0; % for single NLS (test only)

else
    ParamControl.RK4_version = 1; 
    ParamControl.Consider_FWM_or_Not = 1;
    ParamControl.Consider_dFWM_or_Not = 1;
    
end

% PD
ParamControl.Add_PD_noise = 1;
ParamControl.PD_LPF_or_Not = 1;
    ParamControl.PD_case = 3; % 1 PIN 2 APD 3 PIN_TIA 
    
% ADC
ParamControl.Add_ADC_thermal_noise = 1;

%% Rx DSP 

ParamControl.RxDSP_practical_implementation_or_Not =1;

ParamControl.FFT_sharing_or_Not = 1; % FFT sharing in resampling before KK
                                     % and hilbert transform
    ParamControl.real_output_subFFT_num = 1;

ParamControl.Synchronization_or_Not = 1;
ParamControl.Equalization_or_Not = 1;
    ParamControl.PLL_or_Not  = 0;
    ParamControl.Update_Eqtap_or_Not = 1;
    
%% FEC option
ParamControl.FEC_option = 3;% 1 HD-FEC, 2 SD-FEC
                            % 3 SD-FEC 15% OH 2e-2



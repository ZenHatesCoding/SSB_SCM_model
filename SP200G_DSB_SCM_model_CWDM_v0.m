% clc; 
clear;
close all;

%% Initialization
SP200G_DSB_CWDM_Simulator_Control_Param_Init;
SP200G_DSB_CWDM_Simulator_System_Param_Init;

SP_DSB_CWDM_PhysicalLayerSimulator;


disp('BER:')
disp(ParamPhysicalModel.BER_CWDM);

switch ParamControl.FEC_option 
    case 1
        ParamOpt.BER_target = 3.8e-3;
    case 2
        ParamOpt.BER_target = 1.25e-2;
    case 3
        ParamOpt.BER_target = 2e-2;
end

if ParamControl.PCA_Enable_or_Not
    if ParamPhysicalModel.BER_CWDM_avg <= ParamOpt.BER_target
        PCA_Param_Init;

        ParamControl.Plot_Power_Pie_or_Not = 1;
        SP_DSB_Power_Consumption_Analysis;
    else
        disp(['cannot achieve below threshold transmission with ',num2str(ParamRxDSP.LMS_linear_tap),' taps.']);
    end
end
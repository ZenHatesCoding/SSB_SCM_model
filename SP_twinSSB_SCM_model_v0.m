
%% ZP 01/04/2020
% currently only works for single carrier 16 QAM 
clc; 
% clear;
close all;

%% Initialization
SP_twinSSB_Simulator_Control_Param_Init;

SP_twinSSB_Simulator_System_Param_Init;

SP_twinSSB_PhysicalLayerSimulator;

disp('CSPR:')
disp(ParamPhysicalModel.Measured_CSPR);
disp('BER:')
disp(ParamPhysicalModel.BER_avg);

switch ParamControl.FEC_option 
    case 1
        ParamOpt.BER_target = 3.8e-3;
    case 2
        ParamOpt.BER_target = 1.25e-2;
end

if ParamControl.PCA_Enable_or_Not 
    if ParamPhysicalModel.BER_avg <= ParamOpt.BER_target
        PCA_Param_Init;

        ParamControl.Plot_Power_Pie_or_Not = 1;
        SP_SSB_Power_Consumption_Analysis;
    else
        disp(['cannot achieve below threshold transmission with ',num2str(ParamRxDSP.LMS_linear_tap),' taps.']);
    end

    ParamOpt.Vpi_range = 6:6; 
    Psum_Vpi_list = zeros(length(ParamOpt.Vpi_range),1);
    for j1 = 1:length(ParamOpt.Vpi_range)
        ParamSys.Vpi = ParamOpt.Vpi_range(j1);
        SP_SSB_Power_Consumption_Analysis;          
        OptResult.Vpi = ParamSys.Vpi;
        OptResult.Power = ParamPCA.P_sum;
        Psum_Vpi_list(j1,1) = OptResult.Power; 
    end
    Psum_Vpi_list.'
    
    ParamPCA
end
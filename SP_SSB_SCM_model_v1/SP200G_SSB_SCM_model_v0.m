
clc; 
% clear;
close all;

%% Initialization
% 200G means 64 QAM 40 Gbaud
% currently still at 30 Gbaud - under development
SP200G_SSB_Simulator_Control_Param_Init;

SP200G_SSB_Simulator_System_Param_Init;

SP_SSB_PhysicalLayerSimulator;
% while ParamPhysicalModel.BER_avg > 0.45 % let's try again
%     SP_SSB_PhysicalLayerSimulator;
% end
disp('CSPR:')
disp(ParamPhysicalModel.Measured_CSPR);
disp('BER:')
disp(ParamPhysicalModel.BER_avg);
disp('SNR:')
disp(ParamPhysicalModel.SNR_list);
% disp('Tx Loss:')
% disp(ParamPhysicalModel.Tx_optical_Loss_dB);
% disp('ROP:')
% disp(ParamPhysicalModel.Received_Pwr_dBm);

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
        
        
%         ParamOpt.Vpi_range = 1.5; 
%         Psum_Vpi_list = zeros(length(ParamOpt.Vpi_range),1);
%         for j1 = 1:length(ParamOpt.Vpi_range)
%             ParamSys.Vpi = ParamOpt.Vpi_range(j1);
%             SP_SSB_Power_Consumption_Analysis;          
%             OptResult.Vpi = ParamSys.Vpi;
%             OptResult.Power = ParamPCA.P_sum;
%             Psum_Vpi_list(j1,1) = OptResult.Power; 
%         end
%         Psum_Vpi_list.'

        ParamPCA
        a = ParamPCA.P_KK+ParamPCA.P_Rx_resamp+ParamPCA.P_ADC
        b = ParamPCA.P_driver+ParamPCA.P_Tx_resamp+ParamPCA.P_DAC+ParamPCA.P_laser
        c = a+b
    else
        disp(['cannot achieve below threshold transmission with ',num2str(ParamRxDSP.LMS_linear_tap),' taps.']);
    end


end
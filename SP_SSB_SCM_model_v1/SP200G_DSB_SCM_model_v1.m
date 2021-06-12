%% ZP 01/04/2020
% currently only works for single carrier 16 QAM 
clc;
clear;
close all;

%% Initialization
SP200G_DSB_Simulator_Control_Param_Init;

ParamControl.Laser_case = 2;

SP200G_DSB_Simulator_System_Param_Init;
PCA_Param_Init;

%% Optimization
ParamOpt.FFE_tap_range = 3:2:23;
ParamOpt.Vbias_over_Vpi_range = 0.5;
ParamOpt.Vpi_range = 2.5;


switch ParamControl.FEC_option 
    case 1
        ParamOpt.BER_target = 3.8e-3;
    case 2
        ParamOpt.BER_target = 1.25e-2;
    case 3
        ParamOpt.BER_target = 2e-2;
end
ParamOpt.num_trial = 5;

lgd1_idx = 1;
for i4 = 1:length(ParamOpt.Vbias_over_Vpi_range)
    ParamSys.Vbias_over_Vpi = ParamOpt.Vbias_over_Vpi_range(i4);
 
    ParamOpt.Vpp_over_Vpi_range = 0.7:0.05:1;
    ParamOpt.Psum_list = zeros(size(ParamOpt.Vpp_over_Vpi_range));
    ParamOpt.FFE_tap_list = zeros(size(ParamOpt.Vpp_over_Vpi_range));
    
    %*********************************************************************
    for i3 = 1:length(ParamOpt.Vpp_over_Vpi_range)
        ParamSys.Vpp_over_Vpi = ParamOpt.Vpp_over_Vpi_range(i3);
        ParamOpt.Flag_under_FEC_threshold = 0;
        for i2 = 1:length(ParamOpt.FFE_tap_range)    
            ParamRxDSP.LMS_linear_tap = ParamOpt.FFE_tap_range(i2);
            disp(['Vbias/Vpi=',num2str(ParamSys.Vbias_over_Vpi)]);
            disp(['Vpp/Vpi=',num2str(ParamSys.Vpp_over_Vpi)]);
            disp(['numTap=',num2str(ParamRxDSP.LMS_linear_tap)]);
            BER_list = ones(ParamOpt.num_trial,1);
            for i1 = 1:ParamOpt.num_trial
                SP_DSB_PhysicalLayerSimulator;
                idx = 0;
                while ParamPhysicalModel.BER_avg > 0.45 % let's try again
                    idx = idx +1;
                    SP_DSB_PhysicalLayerSimulator;
                end
                BER_list(i1,1) = ParamPhysicalModel.BER_avg;
                
                disp(['trial ',num2str(i1),' finished.']);
            end
            
            ParamOpt.BER_avg_over_Ntrial = mean(BER_list);
        
            if ParamOpt.BER_avg_over_Ntrial <= ParamOpt.BER_target
                ParamOpt.Flag_under_FEC_threshold = 1;
                break;
            end
        end  
        if ParamOpt.Flag_under_FEC_threshold == 0
            disp(['cannot achieve below threshold transmission with ',num2str(ParamOpt.FFE_tap_range(end)),' taps.']);
        else
            disp('Compute Power consumption');
            OptResult.laserpower = ParamLas.laser_power_dBm;
            OptResult.Vbias_over_Vpi = ParamSys.Vbias_over_Vpi;
            OptResult.Vpp_over_Vpi = ParamSys.Vpp_over_Vpi;
            OptResult.num_FFE_tap = ParamRxDSP.LMS_linear_tap;
            
            for j1 = 1:length(ParamOpt.Vpi_range)
                ParamSys.Vpi = ParamOpt.Vpi_range(j1);
                SP_DSB_Power_Consumption_Analysis;          
                OptResult.Vpi = ParamSys.Vpi;
                OptResult.Power = ParamPCA.P_sum;
                ParamOpt.Psum_list(i3) = OptResult.Power; 
                ParamOpt.FFE_tap_list(i3) = OptResult.num_FFE_tap; 
            end
%             figure(1);
%             plot(ParamOpt.Vpi_range,Psum_Vpi_list,'v-','LineWidth',2,'Markersize',12);
%             hold on; grid on;
%             lgd1{lgd1_idx} = [num2str(OptResult.Vbias_over_Vpi),',',num2str(OptResult.Vpp_over_Vpi),...
%                     ',',num2str(OptResult.num_FFE_tap),' taps'];
%             lgd1_idx = lgd1_idx + 1;
            
            if ParamRxDSP.LMS_linear_tap == ParamOpt.FFE_tap_range(1)
                break;
            end
        end
    end
end
% figure(1); legend(lgd1); xlabel('V_\pi [V]'); ylabel('Power Consumption [W]');
% set(gca,'LineWidth',2); set(gca,'FontSize',14);

figure;
subplot(122);
plot(ParamOpt.Vpp_over_Vpi_range,ParamOpt.Psum_list,'x','Markersize',10,'LineWidth',2);
ylabel('Power Consumption [W]');
xlabel('V_{pp}/V_{\pi} [V]'); 

grid on;
set(gca,'LineWidth',2); set(gca,'FontSize',14);


subplot(121);
plot(ParamOpt.Vpp_over_Vpi_range,ParamOpt.FFE_tap_list,'o','Markersize',10,'LineWidth',2);
ylabel('# FFE taps');

xlabel('V_{pp}/V_{\pi}'); 

grid on;
set(gca,'LineWidth',2); set(gca,'FontSize',14);

sgtitle('LaserPwr = 15 dBm, LineWidth = 2 MHz');

disp('optimization finished');
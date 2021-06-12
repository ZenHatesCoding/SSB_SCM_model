% v2, SSFT step convergence test
% clc; 
clear;
close all;

%% Initialization
SP_DSB_CWDM_Simulator_Control_Param_Init;
SP_DSB_CWDM_Simulator_System_Param_Init;

Path = 'D:\xzp\Project\CWDM_NL_model\';
TxT_name = [Path,'13dBmLaser_40km_FullsingleNLS_zeroCD_1310_FFE_study_spacing_1dot5T.txt'];

formatspec = [];
for i = 1:4
    formatspec = [formatspec,'%12.4e '];
    
end
formatspec = [formatspec, '\n'];

for FFE_taps = 11:2:51
    ParamRxDSP.LMS_linear_tap = FFE_taps; 
    
    SP_DSB_CWDM_PhysicalLayerSimulator;
    ParamControl.New_Simulation = 0;
    screen_message = ['FFE_use_',num2str(FFE_taps),'_taps.'];
    disp('BER:')
    disp(ParamPhysicalModel.BER_CWDM);
    
    
    fileID = fopen(TxT_name,'at');

    fprintf(fileID,formatspec,ParamPhysicalModel.BER_CWDM);
    fprintf(fileID,[screen_message,'\n']);
    fclose(fileID);
    disp(screen_message);
end

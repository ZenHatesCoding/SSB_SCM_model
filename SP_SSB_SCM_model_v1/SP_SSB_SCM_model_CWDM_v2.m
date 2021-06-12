% v2, SSFT step convergence test
% clc; 
clear;
close all;

%% Initialization
SP_SSB_CWDM_Simulator_Control_Param_Init;
SP_SSB_CWDM_Simulator_System_Param_Init;

Path = 'D:\xzp\Project\CWDM_NL_model\';
TxT_name = [Path,'16dBmLaser_40km_FullSingleNLS_1550nm_conv_study_spacing_100GHzSpacing.txt'];

formatspec = [];
for i = 1:4
    formatspec = [formatspec,'%12.4e '];
end
formatspec = [formatspec, '\n'];

for numsec = [1 4 10 40 100 200 400 500]
    ParamFib.numsec_SSFT = numsec; 
    SP_SSB_CWDM_PhysicalLayerSimulator;
    screen_message = ['SSFT_use_',num2str(numsec),'_sections.'];
    disp('BER:')
    disp(ParamPhysicalModel.BER_CWDM);
    
    fileID = fopen(TxT_name,'at');

    fprintf(fileID,formatspec,ParamPhysicalModel.BER_CWDM);
    fprintf(fileID,[screen_message,'\n']);
    fclose(fileID);
    disp(screen_message);
end




% v2, SSFT step convergence test
% clc; 
clear;
close all;

%% Initialization
SP_DSB_CWDM_Simulator_Control_Param_Init;
SP_DSB_CWDM_Simulator_System_Param_Init;

Path = 'D:\xzp\Project\CWDM_NL_model\';
TxT_name = [Path,'13dBmLaser_40km_FullSingleNLS_numSEC_400_tap31_zeroCDswp_study_lambda_1296_1305_1315_1324.txt'];

formatspec = [];
for i = 1:4
    formatspec = [formatspec,'%12.4e '];
end
formatspec = [formatspec, '\n'];

%[1.306e-6 1.308e-6 1.31e-6 1.312e-6 1.314e-6];
for lambda_zeroCD = [1.305e-6 1.306e-6 1.308e-6 1.31e-6 1.312e-6 1.314e-6 1.315e-6]
    %% Laser
    ParamLas.Laser_OSNR = 40; % laser OSNR
    ParamLas.Lambda_ref = 1.31e-6; % O band
    ParamLas.omega_ref = 2*pi*ParamGen.c0./ParamLas.Lambda_ref; 
    % ******************************
    % if ParamSig.numChannel == 4
    % %     ParamLas.Lambda0 = [-15.5,-5,5,14.5]*1e-9+ParamLas.Lambda_ref;
    %     ParamLas.Lambda0 = [-15,-5,5,15]*1e-9+ParamLas.Lambda_ref;
    % %     ParamLas.Lambda0 = [-30,-10,10,30]*1e-9+ParamLas.Lambda_ref;
    % else
    %     ParamLas.Lambda0 = ParamLas.Lambda_ref; % O band
    % end
    % 
    % ParamLas.deltaF = ParamGen.c0./ParamLas.Lambda0 - ParamGen.c0./ParamLas.Lambda_ref;
    % ParamLas.omega0 = 2*pi*ParamGen.c0./ParamLas.Lambda0; 
    % ******************************
    ParamLas.deltaF = [2.4,0.8,-0.8,-2.4]*1e12;
    % ParamLas.deltaF = [2.4,1.2,0,-1.2]*1e12;
    % ParamLas.deltaF = [150,50,-50,-150]*1e9;
    ParamLas.omega0 = (ParamLas.deltaF + ParamGen.c0./ParamLas.Lambda_ref) *2*pi;
    ParamLas.Lambda0 = ParamGen.c0./(ParamGen.c0./ParamLas.Lambda_ref + ParamLas.deltaF);

    %% Fiber
    ParamFib.FiberLength = 40e3;
    ParamFib.Lambda_zeroCD = lambda_zeroCD;
    ParamFib.Loss_alpha_dB = 0.35; % [dB/km]
    ParamFib.Loss_alpha = -log(10^(-0.1*ParamFib.Loss_alpha_dB/1e3));
    ParamFib.S0 = 0.089*1e-12/1e-18/1e3; % [ps/nm^2*km]->[s/(m^3)]
    ParamFib.D= ParamFib.S0/4.*(ParamLas.Lambda0-ParamFib.Lambda_zeroCD.^4./ParamLas.Lambda0.^3); 
    ParamFib.D_ref= ParamFib.S0/4.*(ParamLas.Lambda_ref-ParamFib.Lambda_zeroCD.^4./ParamLas.Lambda_ref.^3); 
    ParamFib.ng_ref = 1.4677; 
    ParamFib.ng_ref_lambda = 1.31e-6;

    % D = d Beta1/d lambda; 
    %Beta1 = S0/4(2*pi^2*c^2/w^2 +lam_zeroCD^4*w/4/pi^2/c^2);
    ParamFib.Beta1 = ParamFib.S0/4.*(2*pi^2*ParamGen.c0^2./ParamLas.omega0.^2+...
                     ParamFib.Lambda_zeroCD^4*ParamLas.omega0.^2./(8*pi^2*ParamGen.c0^2));

    ParamFib.Beta1_ng_ref = ParamFib.ng_ref/ParamGen.c0;
    % ParamFib.Beta1_ng_ref = 0;
    ParamFib.omega_ng_ref = 2*pi*ParamGen.c0/ParamFib.ng_ref_lambda;
    ParamFib.Beta1_ng_ref_estimate = ParamFib.S0/4.*(2*pi^2*ParamGen.c0^2./ParamFib.omega_ng_ref.^2+...
                     ParamFib.Lambda_zeroCD^4*ParamFib.omega_ng_ref.^2./(8*pi^2*ParamGen.c0^2));

    ParamFib.Beta1_constant_diff = ParamFib.Beta1_ng_ref-ParamFib.Beta1_ng_ref_estimate;


    ParamFib.Beta1 = ParamFib.Beta1 + ParamFib.Beta1_constant_diff;
    % ParamFib.Beta1 = ParamFib.Beta1 - mean(ParamFib.Beta1);

    ParamFib.Beta1_ref = ParamFib.S0/4.*(2*pi^2*ParamGen.c0^2./ParamLas.omega_ref.^2+...
                     ParamFib.Lambda_zeroCD^4*ParamLas.omega_ref.^2./(8*pi^2*ParamGen.c0^2))+...
                     ParamFib.Beta1_constant_diff;

    % Beta1 = d Beta0/ d w
    ParamFib.Beta0 = ParamFib.S0/4.*(-2*pi^2*ParamGen.c0^2./ParamLas.omega0 +...
        ParamFib.Lambda_zeroCD^4*ParamLas.omega0.^3./(24*pi^2*ParamGen.c0^2))...
        +ParamFib.Beta1_constant_diff.*ParamLas.omega0;

    ParamFib.Beta0 = ParamFib.Beta0 - mean(ParamFib.Beta0);
    ParamFib.Beta2 = -ParamLas.Lambda0.^2.*ParamFib.D./(2*pi*ParamGen.c0);

    ParamFib.Beta2_ref = -ParamLas.Lambda_ref.^2.*ParamFib.D_ref./(2*pi*ParamGen.c0);


    % beta 3 - S0/4(12*pi^2*c^2/w^4+lam_zeroCD^4/4/pi^2/c^2)
    ParamFib.Beta3 = ParamFib.S0/4*(12*pi^2*ParamGen.c0^2./ParamLas.omega0.^4+...
                     ParamFib.Lambda_zeroCD^4/(4*pi^2*ParamGen.c0^2));

    ParamFib.Beta3_ref = ParamFib.S0/4*(12*pi^2*ParamGen.c0^2./ParamLas.omega_ref.^4+...
                     ParamFib.Lambda_zeroCD^4/(4*pi^2*ParamGen.c0^2));              


    ParamFib.Beta4 = ParamFib.S0/4*(-48*pi^2*ParamGen.c0^2./ParamLas.omega0.^5);
    ParamFib.Beta4_ref = ParamFib.S0/4*(-48*pi^2*ParamGen.c0^2./ParamLas.omega_ref.^5);

    ParamFib.Beta5 = ParamFib.S0/4*(240*pi^2*ParamGen.c0^2./ParamLas.omega0.^6);
    ParamFib.Beta5_ref = ParamFib.S0/4*(240*pi^2*ParamGen.c0^2./ParamLas.omega_ref.^6);

    ParamFib.MFD = 9.2e-6; %[m] at 1310 nm              
    ParamFib.Aeff = 1/4*pi*(ParamFib.MFD^2);
    ParamFib.n2 = 2.6e-20; %[m^2/W]
    % ParamFib.gamma = 2*pi*ParamFib.n2/ParamFib.Aeff./ParamLas.Lambda0; % [W^-1/m]
    ParamFib.gamma = 2*pi*ParamFib.n2/ParamFib.Aeff./ParamLas.Lambda_ref; % [W^-1/m]
         % gamma difference between different wavelengths is negligible 
    ParamFib.gamma = ParamFib.gamma *8/9;     
    ParamFib.numsec_SSFT = 400;

    
    %%-------
    SP_DSB_CWDM_PhysicalLayerSimulator;
    screen_message = ['zero_CD_lambda=',num2str(lambda_zeroCD*1e6),'_nm.'];
    disp('BER:')
    disp(ParamPhysicalModel.BER_CWDM);
    
    fileID = fopen(TxT_name,'at');

    fprintf(fileID,formatspec,ParamPhysicalModel.BER_CWDM);
    fprintf(fileID,[screen_message,'\n']);
    fclose(fileID);
    disp(screen_message);
end




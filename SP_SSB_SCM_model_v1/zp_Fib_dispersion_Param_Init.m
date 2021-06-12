ParamFib.Lambda_zeroCD = 1.31e-6; % O band

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
% ParamFib.Beta2_ref = -ParamFib.Beta2_ref;

% beta 3 - S0/4(12*pi^2*c^2/w^4+lam_zeroCD^4/4/pi^2/c^2)
ParamFib.Beta3 = ParamFib.S0/4*(12*pi^2*ParamGen.c0^2./ParamLas.omega0.^4+...
                 ParamFib.Lambda_zeroCD^4/(4*pi^2*ParamGen.c0^2));
             
ParamFib.Beta3_ref = ParamFib.S0/4*(12*pi^2*ParamGen.c0^2./ParamLas.omega_ref.^4+...
                 ParamFib.Lambda_zeroCD^4/(4*pi^2*ParamGen.c0^2));              


ParamFib.Beta4 = ParamFib.S0/4*(-48*pi^2*ParamGen.c0^2./ParamLas.omega0.^5);
ParamFib.Beta4_ref = ParamFib.S0/4*(-48*pi^2*ParamGen.c0^2./ParamLas.omega_ref.^5);

ParamFib.Beta5 = ParamFib.S0/4*(240*pi^2*ParamGen.c0^2./ParamLas.omega0.^6);
ParamFib.Beta5_ref = ParamFib.S0/4*(240*pi^2*ParamGen.c0^2./ParamLas.omega_ref.^6);

ParamFib.MFD = 10.4e-6; %[m] at 1550 nm              
ParamFib.Aeff = 1/4*pi*(ParamFib.MFD^2);
ParamFib.n2 = 2.6e-20; %[m^2/W]

ParamFib.gamma = 2*pi*ParamFib.n2/ParamFib.Aeff./ParamLas.Lambda_ref; % [W^-1/m]
     % gamma difference between different wavelengths is negligible 
ParamFib.gamma = ParamFib.gamma *8/9;  
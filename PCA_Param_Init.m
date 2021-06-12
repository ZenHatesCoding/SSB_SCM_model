% ParamCMOS

% CMOS
ae2 = 0.09365; ae1 = -0.03409; ae0 = 0.02043;
V_7 = 0.6;
EFactor_180 = 26.372;
EFactor_7 = ae2*V_7^2+ae1*V_7+ae0;

% 7 nm scaling from 180 nm
ParamCMOS.E_opG = 402.4*EFactor_7/EFactor_180; %[fJ]
ParamCMOS.E_opR = 2000*EFactor_7/EFactor_180; %[fJ]
ParamCMOS.E_opRO = 7978.2*EFactor_7/EFactor_180; % [fJ]
ParamCMOS.E_opA = 11991*EFactor_7/EFactor_180; %[fJ]
ParamCMOS.E_opM = 95925*EFactor_7/EFactor_180; %[fJ]
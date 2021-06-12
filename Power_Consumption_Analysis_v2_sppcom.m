clear;
PCA_Param_Init;
SP_SSB_Simulator_Control_Param_Init;
SP_SSB_Simulator_System_Param_Init;
ParamControl.FEC_option = 1;
ParamControl.KK_option = 1; 
SP_SSB_PhysicalLayerSimulator;


% Prepare general parameters
switch ParamControl.FEC_option
    case 1
        R_code = 100/112; % FEC + ODU-4
    case 2
        R_code = 100/120;
end
B = ParamSig.Baud_Rate/1e9; %[Gbaud]
MM = 2^ParamSig.SC(1); % 16-QAM
D = ParamFib.D; %[ps/nm-km]
Lkm = ParamFib.FiberLength/1e3; %[km]
n_DAC = ParamDAC.qnbit_DAC;
n_ADC = ParamADC.qnbit_ADC;


%% TxDSP - RC 
n_samp = 2;
n_tap = 32;

nf = 2;
nb = n_tap*log2(MM)/(2*nf*n_samp);
E_opRO_new = nb/n_DAC*ParamCMOS.E_opRO;

ParamECA_DSP.E_RC = 2*n_samp*(n_samp*nf*n_DAC*E_opRO_new+ParamCMOS.E_opA*(nf-1))/(R_code*log2(MM));

%% TxDSP - resample at DAC (upconversion)
n_samp = 64/28;

switch ParamControl.FEC_option
    case 1
        % 896 = 7*128 FFT 
        Nfft1 = 7;
        Nfft2 = 128;
%         N_opCM1 = (Nfft1-1)*(Nfft1-1); 
%         N_opCA1 = Nfft1*(Nfft1-1);
%         N_opM1 = 3*N_opCM1;
%         N_opA1 = 3*N_opCM1 + 2*N_opCA1; % single tap mul 3m + 3a, normal mul 3m + 5a 
        N_opM1 = 8*2;
        N_opA1 = 36*2;
        
        [N_opM2,N_opA2] = numOP_SplitRadix_FFT(Nfft2);

        N_opM = Nfft2*N_opM1 + Nfft1*N_opM2;
        N_opA = Nfft2*N_opA1 + Nfft1*N_opA2;

        % 1024 FFT
        Nfft3 = 1024;
        [N_opM3,N_opA3] = numOP_SplitRadix_FFT(Nfft3);
        N_opM = N_opM + N_opM3;
        N_opA = N_opA + N_opA3;

        % FFT/IFFT, save i/o sequence
        N_opR = 4*(Nfft1*Nfft2+Nfft3);

        
        % ignore gate operation temperarely
        N_opG = 0;
        ParamECA_DSP.E_Tx_resamp = n_samp*(N_opM*ParamCMOS.E_opM+N_opA*ParamCMOS.E_opA+n_DAC*N_opR*ParamCMOS.E_opR+N_opG*ParamCMOS.E_opG*n_DAC)/(Nfft3*R_code*log2(MM));

    case 2
        % 1280 = 5*256 FFT 
        Nfft1 = 5;
        Nfft2 = 256;
%         N_opCM1 = (Nfft1-1)*(Nfft1-1); 
%         N_opCA1 = Nfft1*(Nfft1-1);
%         N_opM1 = 3*N_opCM1;
%         N_opA1 = 3*N_opCM1 + 2*N_opCA1; % single tap mul 3m + 3a, normal mul 3m + 5a 
        N_opM1 = 5*2; 
        N_opA1 = 17*2;
        
        [N_opM2,N_opA2] = numOP_SplitRadix_FFT(Nfft2);

        N_opM = Nfft2*N_opM1 + Nfft1*N_opM2;
        N_opA = Nfft2*N_opA1 + Nfft1*N_opA2;

        % 1024 FFT
        Nfft3 = 1024;
        [N_opM3,N_opA3] = numOP_SplitRadix_FFT(Nfft3);
        N_opM = N_opM + N_opM3;
        N_opA = N_opA + N_opA3;

        % FFT/IFFT, save i/o sequence
        N_opR = 4*(Nfft1*Nfft2+Nfft3);

        
        % ignore gate operation temperarely
        N_opG = 0;
        
        ParamECA_DSP.E_Tx_resamp = n_samp*(N_opM*ParamCMOS.E_opM+N_opA*ParamCMOS.E_opA+n_DAC*N_opR*ParamCMOS.E_opR+N_opG*ParamCMOS.E_opG*n_DAC)/(Nfft1*Nfft2*R_code*log2(MM));

end

%% Tx-DSP Pre-emphasis
n_samp = ParamDAC.DAC_Rate/ParamSig.Baud_Rate;
n_tap = 51;
Nfft = 1024;
Nn = Nfft-n_tap+1;

[N_opM,N_opA] = numOP_SplitRadix_FFT(Nfft);
%single tap mul + FFT + IFFT
N_opM = 2*N_opM + 3*Nfft; 
N_opA = 2*N_opA + 3*Nfft;
N_opR = 2*(Nfft+Nn);

N_opG = 0;
ParamECA_DSP.E_preemp = n_samp*(N_opM*ParamCMOS.E_opM+N_opA*ParamCMOS.E_opA+N_opR*n_DAC*ParamCMOS.E_opR+N_opG*ParamCMOS.E_opG*n_DAC)/(R_code*Nn*log2(MM));

%% Resample after ADC
n_samp = 80/28;
% 1024 point to 1334 point
% 1280 = 256*5  point FFT
Nfft1 = 5;
Nfft2 = 256;

% N_opCM1 = (Nfft1-1)*(Nfft1-1); 
% N_opCA1 = Nfft1*(Nfft1-1);
% N_opM1 = 3*N_opCM1;
% N_opA1 = 3*N_opCM1 + 2*N_opCA1; % single tap mul 3m + 3a, normal mul 3m + 5a 
N_opM1 = 5*2; 
N_opA1 = 17*2;

[N_opM2,N_opA2] = numOP_SplitRadix_FFT(Nfft2);

N_opM = Nfft2*N_opM1 + Nfft1*N_opM2;
N_opA = Nfft2*N_opA1 + Nfft1*N_opA2;

% 1024 point FFT
Nfft4 = 1024;
[N_opM4,N_opA4] = numOP_SplitRadix_FFT(Nfft4);
N_opM = N_opM + N_opM4;
N_opA = N_opA + N_opA4;

% FFT/IFFT, save i/o sequence
N_opR = 4*(Nfft1*Nfft2+Nfft4);
% N_opR = 0;
% ignore gate operation temperarely
N_opG = 0;
ParamECA_DSP.E_Rx_resamp_BeforeKK = n_samp*(N_opM*ParamCMOS.E_opM+N_opA*ParamCMOS.E_opA+n_ADC*N_opR*ParamCMOS.E_opR+N_opG*ParamCMOS.E_opG*n_ADC)/(Nfft1*Nfft2*R_code*log2(MM));

%% Rx-DSP KK
if ParamControl.KK_option == 0
    n_samp = 80/28;
else
    n_samp = 64/28;
end

% hilbert transform
Nfft = 1024;
Nfir = ParamRxDSP.hilbert_tap;

Nn = Nfft-Nfir+1;

[N_opM,N_opA] = numOP_SplitRadix_FFT(Nfft);
% single-tap mul + FFT/IFFT
N_opM = 2*N_opM + 3*Nfft;
N_opA = 2*N_opA + 3*Nfft;

N_opR = 2*(Nfft + Nn)*n_ADC;
N_opG = 0;

%
if ParamControl.KK_option == 0
    N_opM = N_opM +2*Nn; % 2 real mul after cos()/ sin()
    N_opRO = 2*Nfft+2*Nn; % 4 LUT
    N_opG = N_opG + 3*Nn; % 3 switch
else
    N_opM = N_opM + (3*Nn+2*Nfft);
    N_opA = N_opA + Nn*5 + Nn;
    N_opRO = Nfft; % 1 LUT
    N_opG = N_opG + 2*Nn + Nfft; % 2 switch
    N_opR = 2*N_opR + 3*2*n_ADC*Nfft; % divide/multiplied by 2 
end

ParamECA_DSP.E_KK =  n_samp*(N_opM*ParamCMOS.E_opM+N_opA*ParamCMOS.E_opA+...
    N_opRO*n_ADC*ParamCMOS.E_opRO+N_opR*ParamCMOS.E_opR+N_opG*ParamCMOS.E_opG*n_ADC)/(R_code*Nn*log2(MM));

% ParamECA_DSP.E_KK = n_samp*((Nfir/2+2)*ParamCMOS.E_opM+(Nfir/2-1)*ParamCMOS.E_opA)/(R_code*log2(MM));
%% Rx-DSP Resample at Rx (downconversion + LPF)
if ParamControl.KK_option == 0
    n_samp = 2;
    % 1536 = 3*512  point FFT
    Nfft1 = 5;
    Nfft2 = 256;

%     N_opCM1 = (Nfft1-1)*(Nfft1-1); 
%     N_opCA1 = Nfft1*(Nfft1-1);
%     N_opM1 = 3*N_opCM1;
%     N_opA1 = 3*N_opCM1 + 2*N_opCA1; % single tap mul 3m + 3a, normal mul 3m + 5a 
    N_opM1 = 5*2;
    N_opA1 = 17*2;
    
    [N_opM2,N_opA2] = numOP_SplitRadix_FFT(Nfft2);

    N_opM = Nfft2*N_opM1 + Nfft1*N_opM2;
    N_opA = Nfft2*N_opA1 + Nfft1*N_opA2;
    
    Nfft3 = 7;
    Nfft4 = 128;
%     N_opCM3 = (Nfft3-1)*(Nfft3-1); 
%     N_opCA3 = Nfft3*(Nfft3-1);
%     N_opM3 = 3*N_opCM3;
%     N_opA3 = 3*N_opCM3 + 2*N_opCA3; % single tap mul 3m + 3a, normal mul 3m + 5a 
    N_opM3 = 8*2;
    N_opA3 = 36*2;

    [N_opM4,N_opA4] = numOP_SplitRadix_FFT(Nfft4);

    N_opM = N_opM + Nfft3*N_opM4 + Nfft4*N_opM3;
    N_opA = N_opA + Nfft3*N_opA4 + Nfft4*N_opA3;
% %     % 1024 point FFT
% %     Nfft3 = 1024;
% %     [N_opM3,N_opA3] = numOP_SplitRadix_FFT(Nfft3);
% %     N_opM = N_opM + N_opM3;
% %     N_opA = N_opA + N_opA3;

    % FFT/IFFT, save i/o sequence
    N_opR = 4*(Nfft1*Nfft2+Nfft3*Nfft4);


    % ignore gate operation temperarely
    N_opG = 0;
    ParamECA_DSP.E_Rx_resamp = n_samp*(N_opM*ParamCMOS.E_opM+N_opA*ParamCMOS.E_opA+n_ADC*N_opR*ParamCMOS.E_opR+N_opG*ParamCMOS.E_opG*n_ADC)/(Nfft3*Nfft4*R_code*log2(MM));
else
    ParamECA_DSP.E_Rx_resamp = ParamECA_DSP.E_Tx_resamp;
end
%% Rx-DSP CDC
n_samp = 2;
Ncd = round(0.032*B^2*Lkm*D/1000);
Nfft = 1024;
Nn = Nfft-Ncd+1;

[N_opM,N_opA] = numOP_SplitRadix_FFT(Nfft);
% single-tap mul + FFT/IFFT
N_opM = 2*N_opM + 3*Nfft;
N_opA = 2*N_opA + 3*Nfft;

N_opR = 2*(Nfft + Nn);
N_opG = 0;

nb = 10; 
E_opM_new = ParamCMOS.E_opM/(n_ADC^2)*(nb^2);
E_opA_new = ParamCMOS.E_opA/n_ADC*nb;

ParamECA_DSP.E_CDC =  n_samp*(N_opM*E_opM_new+N_opA*E_opA_new+...
    N_opRO*n_ADC*ParamCMOS.E_opRO+N_opR*n_ADC*ParamCMOS.E_opR+N_opG*ParamCMOS.E_opG*n_ADC)/(R_code*Nn*log2(MM));

if ParamControl.CD_Compensation_or_Not == 0
    ParamECA_DSP.E_CDC = 0;
end
%% RxDSP - Time Recovery
n_samp = 2;
B_TR = 256;
N_opM = 7*B_TR+2;
N_opA = 28*B_TR+1;
N_opRO = 4.5*n_ADC*B_TR;
N_opR = n_ADC*(16*B_TR+6);
N_opG = 5*n_ADC*B_TR;

ParamECA_DSP.E_TR = n_samp*(N_opM*ParamCMOS.E_opM+N_opA*ParamCMOS.E_opA+N_opRO*ParamCMOS.E_opRO+N_opR*ParamCMOS.E_opR...
       +N_opG*ParamCMOS.E_opG)/(R_code*B_TR*log2(MM));
   
%% FFE -Time Domain
n_samp = 1;
n_tap = ParamRxDSP.LMS_linear_tap;
N_opM = 3*(2*n_tap+1);
N_opA = 5*(2*n_tap) + 3 + 2*2*(n_tap-1) + 2 + 2;
N_opR = 2*(n_tap + n_tap + 1);
N_opG = 0;
ParamECA_DSP.E_FFE_adpt = n_samp*(N_opM*ParamCMOS.E_opM+N_opA*ParamCMOS.E_opA+N_opR*ParamCMOS.E_opR...
       +N_opG*ParamCMOS.E_opG)/(R_code*log2(MM));

%% Convert DSP Ener-consump to Pwr-consump
ParamSig.Info_Bit_Rate = 100e9; 
ParamPCA.P_RC = ParamECA_DSP.E_RC*1e-15*ParamSig.Info_Bit_Rate;
ParamPCA.P_Tx_resamp = ParamECA_DSP.E_Tx_resamp*1e-15*ParamSig.Info_Bit_Rate;
ParamPCA.P_preemp = ParamECA_DSP.E_preemp*1e-15*ParamSig.Info_Bit_Rate;
ParamPCA.P_KK = ParamECA_DSP.E_KK*1e-15*ParamSig.Info_Bit_Rate;
ParamPCA.P_Rx_resamp_BeforeKK = ParamECA_DSP.E_Rx_resamp_BeforeKK*1e-15*ParamSig.Info_Bit_Rate;
ParamPCA.P_Rx_resamp = ParamECA_DSP.E_Rx_resamp*1e-15*ParamSig.Info_Bit_Rate;
ParamPCA.P_CDC = ParamECA_DSP.E_CDC*1e-15*ParamSig.Info_Bit_Rate;
ParamPCA.P_TR = ParamECA_DSP.E_TR*1e-15*ParamSig.Info_Bit_Rate;
ParamPCA.P_FFE_adpt = ParamECA_DSP.E_FFE_adpt*1e-15*ParamSig.Info_Bit_Rate;

%% Hardware
if ParamControl.VSB_or_Not
    num_DAC = 1;
    num_Mod = 1;
    num_driver = 1;
else
    num_DAC = 2;
    num_Mod = 2;
    num_driver = 2;
end
num_ADC = 1;

ParamPCA.P_DAC = num_DAC*2.7e-12/80e9*ParamDAC.DAC_Rate/7*ParamDAC.qnbit_DAC*ParamSig.Info_Bit_Rate ;
% ParamPCA.P_ADC = num_ADC*5e-12/80e9*ParamADC.ADC_Rate/7*ParamADC.qnbit_ADC*ParamSig.Info_Bit_Rate ;
ParamPCA.P_ADC = 0.303/60/7*8*(84-64);
if ParamControl.FEC_option == 1
    ParamPCA.P_FEC = (0.42*0.05+0.42*0.95*(ParamSig.Info_Bit_Rate/400e9))/2;
else
    ParamPCA.P_FEC = 0.42/4;
end
ParamPCA.P_PD = 0.169 + 10^(ParamPhysicalModel.Received_Pwr_dB/10)/1e3*0.76*5;

ParamPhysicalModel.ModVpp = ParamSys.Vpi*ParamSys.Vpp_over_Vpi;
ParamPhysicalModel.ModVrms = mean([ParamPhysicalModel.Vrms_over_VM_I,ParamPhysicalModel.Vrms_over_VM_Q])*ParamPhysicalModel.ModVpp/2;

ParamPhysicalModel.ModVbias = ParamSys.Vpi*ParamSys.Vbias_over_Vpi;
ParamPCA.P_mod = num_Mod*((ParamPhysicalModel.ModVrms)^2/50 + (ParamPhysicalModel.ModVbias)^2/50 );
% ParamPCA.P_mod = num_Mod*(ParamPhysicalModel.ModVrms)^2/50;


ParamPCA.P_driver = num_driver*0.54/3*ParamPhysicalModel.ModVpp;

ParamPCA.P_laser_TEC = 1.2*2.4;
switch ParamControl.Laser_case
    case 1
        ParamPCA.P_laser = 0.08*1.8;
    case 2
        ParamPCA.P_laser = 0.15*1.8;
    case 3
        ParamPCA.P_laser = 0.3*2.5;
end


P_all = [ParamPCA.P_DAC,ParamPCA.P_ADC,...
         ParamPCA.P_RC,ParamPCA.P_Tx_resamp,...
         ParamPCA.P_preemp, ParamPCA.P_KK,...
         ParamPCA.P_Rx_resamp, ParamPCA.P_CDC,...
         ParamPCA.P_TR, ParamPCA.P_FFE_adpt,...
         ParamPCA.P_FEC,...
         ParamPCA.P_laser_TEC,ParamPCA.P_laser...
         ParamPCA.P_mod,ParamPCA.P_driver,...
         ParamPCA.P_PD];
ParamPCA.P_sum = sum(P_all);
% disp(ParamPCA);

if ParamControl.Plot_Power_Pie_or_Not
    label_all = {'DAC','ADC','RC','resampTx','pre-emp','KK',...
                 'resampRx','CDC','TR','FFE-adapt',...
                 'FEC','laser(TEC)','laser(tunable)','modulator','driver','PIN-TIA'};
    explode = ones(size(P_all));
    figure; pie(P_all,explode,label_all);

end
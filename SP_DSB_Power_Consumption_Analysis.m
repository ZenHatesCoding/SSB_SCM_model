% Prepare general parameters
switch ParamControl.FEC_option
    case 1
        R_code = 100/112; % FEC + ODU-4
    case 2
        R_code = 100/120;
    case 3
        R_code = 100/115;
end

MM = 2^ParamSig.SC(1); 

n_DAC = ParamDAC.qnbit_DAC;
n_ADC = ParamADC.qnbit_ADC;


%% TxDSP - RC 
n_samp = 2;
n_tap = 32;

nf = 2;
nb = n_tap*ceil(log2(MM))/(nf*n_samp);
E_opRO_new = nb/n_DAC*ParamCMOS.E_opRO;

ParamECA_DSP.E_RC = n_samp*(n_samp*nf*n_DAC*E_opRO_new+ParamCMOS.E_opA*(nf-1))/(R_code*log2(MM));

%% TxDSP - resample at DAC (upconversion)
n_samp = ParamDAC.DAC_Rate/ParamSig.Baud_Rate;

switch ParamControl.FEC_option
    case 1
        % 896 = 7*128 FFT 
        Nfft1 = 7;
        Nfft2 = 256;
        
        N_opM1 = 8*2;
        N_opA1 = 36*2;
        
        [N_opM2,N_opA2] = numOP_SplitRadix_FFT(Nfft2);

        N_opM = Nfft2*N_opM1 + Nfft1*N_opM2;
        N_opA = Nfft2*N_opA1 + Nfft1*N_opA2;
        if ParamControl.FFT_sharing_or_Not == 1
            N_opM = N_opM/2;
            N_opA = N_opA/2 + Nfft1*Nfft2*4/2;
        end
        % 1024 FFT
        Nfft3 = 1024;
        [N_opM3,N_opA3] = numOP_SplitRadix_FFT(Nfft3);
        % n-fft of real sequence can be implemented as n/2-fft of complex
        % sequence
        if ParamControl.FFT_sharing_or_Not == 1
            N_opM3 = N_opM3/2;
            N_opA3 = N_opA3/2 + Nfft3*4/2;
        end
        
        N_opM = N_opM + N_opM3;
        N_opA = N_opA + N_opA3;

        % FFT/IFFT, save i/o sequence
        if ParamControl.FFT_sharing_or_Not == 0
            N_opR = 3*(Nfft1*Nfft2+Nfft3);
        else
            N_opR = 4*(Nfft1*Nfft2+Nfft3)/2;
        end
        % ignore gate operation temperarely
        N_opG = 0;
        ParamECA_DSP.E_Tx_resamp = n_samp*(N_opM*ParamCMOS.E_opM+N_opA*ParamCMOS.E_opA+n_DAC*N_opR*ParamCMOS.E_opR+N_opG*ParamCMOS.E_opG*n_DAC)/(Nfft3*R_code*log2(MM));

    case 2
        % 1280 = 5*256 FFT 
        Nfft1 = 5;
        Nfft2 = 128;

        N_opM1 = 5*2; 
        N_opA1 = 17*2;
        
        [N_opM2,N_opA2] = numOP_SplitRadix_FFT(Nfft2); 
        % n-fft of real sequence can be implemented as n/2-fft of complex
        % sequence
        

        N_opM = Nfft2/2*N_opM1 + Nfft1*N_opM2;
        N_opA = Nfft2/2*N_opA1 + Nfft1*N_opA2;
        
        if ParamControl.FFT_sharing_or_Not == 1
            N_opM = N_opM/2;
            N_opA = N_opA/2 + Nfft1*Nfft2*4/2;
        end

        % 1024 FFT
        Nfft3 = 1024;
        [N_opM3,N_opA3] = numOP_SplitRadix_FFT(Nfft3);
        if ParamControl.FFT_sharing_or_Not == 1
            N_opM3 = N_opM3/2;
            N_opA3 = N_opA3/2 + Nfft3*4/2;
        end
        N_opM = N_opM + N_opM3;
        N_opA = N_opA + N_opA3;


        % FFT/IFFT, save i/o sequence
        if ParamControl.FFT_sharing_or_Not == 0
            N_opR = 3*(Nfft1*Nfft2+Nfft3);
        else
            N_opR = 4*(Nfft1*Nfft2+Nfft3)/2;
        end

        
        % ignore gate operation temperarely
        N_opG = 0;
        
        ParamECA_DSP.E_Tx_resamp = n_samp*(N_opM*ParamCMOS.E_opM+N_opA*ParamCMOS.E_opA+n_DAC*N_opR*ParamCMOS.E_opR+N_opG*ParamCMOS.E_opG*n_DAC)/(Nfft1*Nfft2*R_code*log2(MM));
    case 3
        % 640 = 3*256 FFT 
        Nfft1 = 3;
        Nfft2 = 256;
        
        N_opM1 = 2*2; 
        N_opA1 = 6*2;
        
        [N_opM2,N_opA2] = numOP_SplitRadix_FFT(Nfft2);

        N_opM = Nfft2*N_opM1 + Nfft1*N_opM2;
        N_opA = Nfft2*N_opA1 + Nfft1*N_opA2;
        if ParamControl.FFT_sharing_or_Not == 1
            N_opM = N_opM/2;
            N_opA = N_opA/2 + Nfft1*Nfft2*4/2;
        end
        % 1024 FFT
        Nfft3 = 1024;
        [N_opM3,N_opA3] = numOP_SplitRadix_FFT(Nfft3);
        % n-fft of real sequence can be implemented as n/2-fft of complex
        % sequence
        if ParamControl.FFT_sharing_or_Not == 1
            N_opM3 = N_opM3/2;
            N_opA3 = N_opA3/2 + Nfft3*4/2;
        end
        
        N_opM = N_opM + N_opM3;
        N_opA = N_opA + N_opA3;

        % FFT/IFFT, save i/o sequence
        if ParamControl.FFT_sharing_or_Not == 0
            N_opR = 3*(Nfft1*Nfft2+Nfft3);
        else
            N_opR = 4*(Nfft1*Nfft2+Nfft3)/2;
        end
        % ignore gate operation temperarely
        N_opG = 0;
        ParamECA_DSP.E_Tx_resamp = n_samp*(N_opM*ParamCMOS.E_opM+N_opA*ParamCMOS.E_opA+n_DAC*N_opR*ParamCMOS.E_opR+N_opG*ParamCMOS.E_opG*n_DAC)/(Nfft1*Nfft2*R_code*log2(MM));
 
end

%% Tx-DSP Pre-emphasis 
n_samp = ParamDAC.DAC_Rate/ParamSig.Baud_Rate;
% n_tap = floor(51*ParamDAC.DAC_Rate/64e9);
n_tap = 32;
Nfft = 1024;
Nn = Nfft-n_tap+1;

[N_opM,N_opA] = numOP_SplitRadix_FFT(Nfft);
if ParamControl.FFT_sharing_or_Not == 1
    N_opM = N_opM/2;
    N_opA = N_opA/2 + Nfft*4/2;
end 
%single tap mul + FFT + IFFT
N_opM = 2*N_opM + 3*Nfft; 
N_opA = 2*N_opA + 3*Nfft;
if ParamControl.FFT_sharing_or_Not == 0
    N_opR = 3*(Nfft+Nn);
else
    N_opR = 4*(Nfft+Nn)/2; 
end


N_opG = 0;
ParamECA_DSP.E_preemp = n_samp*(N_opM*ParamCMOS.E_opM+N_opA*ParamCMOS.E_opA+N_opR*n_DAC*ParamCMOS.E_opR+N_opG*ParamCMOS.E_opG*n_DAC)/(R_code*Nn*log2(MM));



%% Rx-DSP Resample at Rx 
if ParamADC.ADC_Rate == 64e9
    n_samp = 2;
    Nfft1 = 1024;

    [N_opM1,N_opA1] = numOP_SplitRadix_FFT(Nfft1);
    if ParamControl.FFT_sharing_or_Not == 1
        N_opM1 = N_opM1/2;
        N_opA1 = N_opA1/2 + Nfft1*4/2;
    end

    Nfft2 = 7;
    Nfft3 = 256;

    N_opM2 = 8*2;
    N_opA2 = 36*2;

    [N_opM3,N_opA3] = numOP_SplitRadix_FFT(Nfft3);
    N_opM = Nfft2*N_opM3 + Nfft3*N_opM2;
    N_opA = Nfft2*N_opA3 + Nfft3*N_opA2;
    if ParamControl.FFT_sharing_or_Not == 1
        N_opM = N_opM/2;
        N_opA = N_opA/2 + Nfft2*Nfft3*4/2;
    end
    N_opM = N_opM1 + N_opM;
    N_opA = N_opA1 + N_opA;

    % FFT/IFFT, save i/o sequence
    if ParamControl.FFT_sharing_or_Not == 0
        N_opR = 3*Nfft1 + 3*Nfft2*Nfft3;
    else
        N_opR = 4*(Nfft1 + Nfft2*Nfft3)/2;
    end
    % ignore gate operation temperarely
    N_opG = 0;
    ParamECA_DSP.E_Rx_resamp = n_samp*(N_opM*ParamCMOS.E_opM+N_opA*ParamCMOS.E_opA+n_ADC*N_opR*ParamCMOS.E_opR+N_opG*ParamCMOS.E_opG*n_ADC)/(Nfft2*Nfft3*R_code*log2(MM));
elseif ParamADC.ADC_Rate == 75e9
    %1280 - 1024
    n_samp = 2;
    Nfft1 = 5;
    Nfft2 = 128;
    
    N_opM1 = 5*2; 
    N_opA1 = 17*2;
    [N_opM2,N_opA2] = numOP_SplitRadix_FFT(Nfft2);

    N_opM = Nfft2*N_opM1 + Nfft1*N_opM2;
    N_opA = Nfft2*N_opA1 + Nfft1*N_opA2;
    if ParamControl.FFT_sharing_or_Not == 1
        N_opM = N_opM/2;
        N_opA = N_opA/2 + Nfft1*Nfft2*4/2;
    end


    Nfft3 = 1024;
    [N_opM3,N_opA3] = numOP_SplitRadix_FFT(Nfft3);
    if ParamControl.FFT_sharing_or_Not == 1
        N_opM3 = N_opM3/2;
        N_opA3 = N_opA3/2 + Nfft3*4/2;
    end
    N_opM = N_opM + N_opM3;
    N_opA = N_opA + N_opA3;

    % FFT/IFFT, save i/o sequence
    if ParamControl.FFT_sharing_or_Not == 0
        N_opR = 3*Nfft1*Nfft2 + 3*Nfft3;
    else
        N_opR = 4*(Nfft1*Nfft2 + Nfft3)/2;
    end
    % ignore gate operation temperarely
    N_opG = 0;
    ParamECA_DSP.E_Rx_resamp = n_samp*(N_opM*ParamCMOS.E_opM+N_opA*ParamCMOS.E_opA+n_ADC*N_opR*ParamCMOS.E_opR+N_opG*ParamCMOS.E_opG*n_ADC)/(Nfft3*R_code*log2(MM));
elseif ParamADC.ADC_Rate == 138e9
    n_samp = 2;
    Nfft1 = 1024;

    [N_opM1,N_opA1] = numOP_SplitRadix_FFT(Nfft1);
    if ParamControl.FFT_sharing_or_Not == 1
        N_opM1 = N_opM1/2;
        N_opA1 = N_opA1/2 + Nfft1*4/2;
    end

    Nfft2 = 3;
    Nfft3 = 256;

    N_opM2 = 2*2; 
    N_opA2 = 6*2;

    [N_opM3,N_opA3] = numOP_SplitRadix_FFT(Nfft3);
    N_opM = Nfft2*N_opM3 + Nfft3*N_opM2;
    N_opA = Nfft2*N_opA3 + Nfft3*N_opA2;
    if ParamControl.FFT_sharing_or_Not == 1
        N_opM = N_opM/2;
        N_opA = N_opA/2 + Nfft2*Nfft3*4/2;
    end
    N_opM = N_opM1 + N_opM;
    N_opA = N_opA1 + N_opA;

    % FFT/IFFT, save i/o sequence
    if ParamControl.FFT_sharing_or_Not == 0
        N_opR = 3*Nfft1 + 3*Nfft2*Nfft3;
    else
        N_opR = 4*(Nfft1 + Nfft2*Nfft3)/2;
    end
    % ignore gate operation temperarely
    N_opG = 0;
    ParamECA_DSP.E_Rx_resamp = n_samp*(N_opM*ParamCMOS.E_opM+N_opA*ParamCMOS.E_opA+n_ADC*N_opR*ParamCMOS.E_opR+N_opG*ParamCMOS.E_opG*n_ADC)/(Nfft1*R_code*log2(MM));

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

% ParamECA_DSP.E_TR = n_samp*(N_opM*ParamCMOS.E_opM+N_opA*ParamCMOS.E_opA+N_opRO*ParamCMOS.E_opRO+N_opR*ParamCMOS.E_opR...
%        +N_opG*ParamCMOS.E_opG)/(R_code*B_TR*4);  
%% FFE -Time Domain
n_samp = 1;
% n_tap = ParamRxDSP.LMS_linear_tap;
n_tap = 31;
N_opM = 2*n_tap;
N_opA = n_tap + (n_tap-1) + 1; % n_tap-1: FFE, n_tap: update, 1:error compute
N_opR = 2*(n_tap + n_tap);
N_opG = 0;
ParamECA_DSP.E_FFE_adpt = n_samp*(N_opM*ParamCMOS.E_opM+N_opA*ParamCMOS.E_opA+N_opR*ParamCMOS.E_opR...
       +N_opG*ParamCMOS.E_opG)/(R_code*log2(MM));

%% Convert DSP Ener-consump to Pwr-consump
ParamSig.Info_Bit_Rate = 100e9;
if ParamControl.FEC_option == 3
    ParamSig.Info_Bit_Rate = 200e9;
end
ParamPCA.P_RC = ParamECA_DSP.E_RC*1e-15*ParamSig.Info_Bit_Rate;
ParamPCA.P_Tx_resamp = ParamECA_DSP.E_Tx_resamp*1e-15*ParamSig.Info_Bit_Rate; % complex
ParamPCA.P_preemp = ParamECA_DSP.E_preemp*1e-15*ParamSig.Info_Bit_Rate; %complex
ParamPCA.P_Rx_resamp_beforeKK = 0;
ParamPCA.P_KK = 0;
ParamPCA.P_Rx_resamp = ParamECA_DSP.E_Rx_resamp*1e-15*ParamSig.Info_Bit_Rate;
ParamPCA.P_CDC = 0;
ParamPCA.P_TR = ParamECA_DSP.E_TR*1e-15*ParamSig.Info_Bit_Rate;
ParamPCA.P_FFE_adpt = ParamECA_DSP.E_FFE_adpt*1e-15*ParamSig.Info_Bit_Rate;
%% Hardware
if ParamControl.IM_or_Not
    num_DAC = 1;
    num_Mod = 1;
    num_driver = 1;
else
    num_DAC = 2;
    num_Mod = 2;
    num_driver = 2;
end
num_ADC = 1;

ParamPCA.P_DAC = num_DAC*2.7e-1/80e9*ParamDAC.DAC_Rate/7*ParamDAC.qnbit_DAC ;
ParamPCA.P_ADC = num_ADC*5e-1/80e9*ParamADC.ADC_Rate/7*ParamADC.qnbit_ADC ;

if ParamOpt.BER_target == 3.8e-3
    ParamPCA.P_FEC = (0.42*0.05+0.42*0.95*(ParamSig.Info_Bit_Rate/400e9))/2;
else
    ParamPCA.P_FEC = 0.42/4;
end

if ParamControl.FEC_option == 3 % 200 G
    ParamPCA.P_FEC = ParamPCA.P_FEC*2;
end

% ParamPhysicalModel.Received_Pwr_dBm = -3;
switch ParamControl.PD_case 
    case 1 
        ParamPCA.P_PD = 10^(ParamPhysicalModel.Received_Pwr_dBm/10)/1e3*ParamPD.R_pin *ParamPD.Vbias;
    case 3
        ParamPCA.P_PD = 0.169 + 10^(ParamPhysicalModel.Received_Pwr_dBm/10)/1e3*ParamPD.R_pin *ParamPD.Vbias;
end

ParamPhysicalModel.ModVpp = ParamSys.Vpi*ParamSys.Vpp_over_Vpi;
if ParamControl.IM_or_Not
    ParamPhysicalModel.ModVrms = ParamPhysicalModel.Vrms_over_VM*ParamPhysicalModel.ModVpp/2;
else
    ParamPhysicalModel.ModVrms = mean([ParamPhysicalModel.Vrms_over_VM_I,ParamPhysicalModel.Vrms_over_VM_Q])*ParamPhysicalModel.ModVpp/2;
end
ParamPhysicalModel.ModVbias = ParamSys.Vpi*ParamSys.Vbias_over_Vpi;
ParamPCA.P_mod = num_Mod*((ParamPhysicalModel.ModVrms)^2/50 + (ParamPhysicalModel.ModVbias)^2/50 );
% ParamPCA.P_mod = num_Mod*(ParamPhysicalModel.ModVrms)^2/50;


ParamPCA.P_driver = num_driver*0.54/3*ParamPhysicalModel.ModVpp;
if ParamControl.TEC_or_Not
    ParamPCA.P_laser_TEC = 1.2*2.4;
else
    ParamPCA.P_laser_TEC = 0;
end


if ParamLas.laser_power_dBm == 10
    ParamPCA.P_laser = 0.08*1.8;
elseif ParamLas.laser_power_dBm == 13
    ParamPCA.P_laser = 0.15*1.8;
elseif ParamLas.laser_power_dBm == 16
    ParamPCA.P_laser = 0.3*2.5;
elseif ParamLas.laser_power_dBm < 16 && ParamLas.laser_power_dBm > 13
    ParamLas.laser_power = 10^(0.1*ParamLas.laser_power_dBm)/1e3;
    ParamPCA.P_laser = 0.15*1.8 + (ParamLas.laser_power - 0.02)/(0.04-0.02)*(0.3*2.5-0.15*1.8); 
end


     
P_all = [ParamPCA.P_DAC,ParamPCA.P_ADC,...
         ParamPCA.P_RC,ParamPCA.P_Tx_resamp,...
         ParamPCA.P_preemp, ...
         ParamPCA.P_Rx_resamp_beforeKK, ParamPCA.P_KK,...
         ParamPCA.P_Rx_resamp, ParamPCA.P_CDC,...
         ParamPCA.P_TR, ParamPCA.P_FFE_adpt,...
         ParamPCA.P_FEC,...
         ParamPCA.P_laser_TEC,ParamPCA.P_laser...
         ParamPCA.P_mod,ParamPCA.P_driver,...
         ParamPCA.P_PD];
ParamPCA.P_sum = sum(P_all);
% disp(ParamPCA);

if ParamControl.Plot_Power_Pie_or_Not
    label_all = {'DAC','ADC','RC','resampTx','pre-emp','resampbeforeKK','KK',...
                 'resampRx','CDC','TR','FFE-adapt',...
                 'FEC','laser(TEC)','laser(tunable)','modulator','driver','PD'};
    explode = ones(size(P_all));
    figure; pie(P_all,explode,label_all);
    

end



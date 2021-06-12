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

f_clk  = 2e9;

%% TxDSP - Modulation 
n_samp = 2;
n_tap = 1;

nf = 2;
nb = n_tap*log2(MM)/(2*nf*n_samp);
E_opRO_new = nb/n_DAC*ParamCMOS.E_opRO;

ParamECA_DSP.E_RC = 2*n_samp*(n_samp*nf*n_DAC*E_opRO_new+ParamCMOS.E_opA*(nf-1))/(R_code*log2(MM));

%% TxDSP - resample at DAC + upconversion + Pre-emphasis
n_samp = ParamDAC.DAC_Rate/ParamSig.Baud_Rate;
N_h = 14*ParamDAC.DAC_Rate/ParamSig.Baud_Rate/2; %overlap-and-save
switch ParamControl.FEC_option
    case 1
        % 896 = 7*128 FFT 
        Nfft1 = 7;
        Nfft2 = 128;
        
        N_opM1 = 8*2;
        N_opA1 = 36*2;
        
        [N_opM2,N_opA2] = numOP_SplitRadix_FFT(Nfft2);

        N_opM = Nfft2*N_opM1 + Nfft1*N_opM2;
        N_opA = Nfft2*N_opA1 + Nfft1*N_opA2;

        % 1024 FFT
        Nfft3 = 1024;
        [N_opM3,N_opA3] = numOP_SplitRadix_FFT(Nfft3);
        % n-fft of real sequence can be implemented as n/2-fft of complex
        % sequence
        
        if ParamControl.VSB_or_Not == 1
            N_opM3 = N_opM3/2; 
            N_opA3 = N_opA3/2 + Nfft3/2;
            N_opR = 4*Nfft1*Nfft2+3*(Nfft3-N_h);
        else
            N_opA3 = N_opA3 + Nfft3;
            N_opR = 6*Nfft1*Nfft2+6*(Nfft3-N_h);
        end
        
        N_tap = ceil(Nfft3/ParamDAC.DAC_Rate*ParamSig.Baud_Rate*(1+ParamSig.roll_off));
        
        N_opM = N_opM + N_opM3 + 3*N_tap;
        N_opA = N_opA + N_opA3 + 3*N_tap;

        % FFT/IFFT, save i/o sequence
        
        
        
        
        n_parallel = ParamDAC.DAC_Rate/f_clk;
        
        if ParamControl.VSB_or_Not == 1
            N_opG = (2*Nfft1*Nfft2+2*(Nfft3-N_h))*ceil(log2(n_parallel));
        else
            N_opG = 2*(Nfft1*Nfft2+2*(Nfft3-N_h))*ceil(log2(n_parallel));
        end
        
        
        ParamECA_DSP.E_Tx_resamp = n_samp*(N_opM*ParamCMOS.E_opM+N_opA*ParamCMOS.E_opA+n_DAC*N_opR*ParamCMOS.E_opR+N_opG*ParamCMOS.E_opG*n_DAC)/...
            ((Nfft3-N_h)*R_code*log2(MM));

end



%% Rx-DSP resampling before KK
% adjust ADC sampling rate directly


%% Rx-DSP KK
n_samp = ParamRxDSP.KKoverSamp;

% hilbert transform
Nfft = 1024;

[N_opM,N_opA] = numOP_SplitRadix_FFT(Nfft);
% single-tap mul + FFT/IFFT
% multiplied by ii/-1i, swtich real and imag
N_opM = 2*N_opM;
N_opA = 2*N_opA;


n_parallel = ParamRxDSP.KKoverSamp*ParamSig.Baud_Rate/f_clk;
N_opM = N_opM/2;
N_opA = N_opA/2 + 3*Nfft;

if ParamControl.KK_option == 3
    N_opR = (3*Nfft+Nfft-ParamRxDSP.hilbert_tap)*n_ADC + Nfft;
else
    N_opR = (3*Nfft+Nfft-ParamRxDSP.hilbert_tap+2*Nfft)*n_ADC+Nfft;
end
N_opG = (Nfft+2*(Nfft-ParamRxDSP.hilbert_tap))*ceil(log2(n_parallel));

% up to now Hilbert transform considered

if ParamControl.KK_option == 0       
    N_opM = N_opM+ Nfft; % 1 real mul after cos()
    N_opRO = 3*Nfft; % 3 LUT
   
    

    
elseif  ParamControl.KK_option == 1
    N_opM = N_opM + 3*Nfft;
    N_opA = N_opA + 2*Nfft;
    N_opR = N_opR + 3*2*Nfft*n_ADC + Nfft*2;
    N_opRO = Nfft; % 1 LUT

elseif  ParamControl.KK_option == 2 % share FFT
    N_opM = N_opM + Nfft;
    N_opA = N_opA + Nfft;
    N_opRO = 2*Nfft; 
    N_opR = N_opR + Nfft*2; % bit flipping
elseif  ParamControl.KK_option == 3 % share FFT % don't use Hilbert transform
    N_opM = N_opM + Nfft + 3*Nfft;
    N_opA = N_opA + Nfft + 3*Nfft;
    N_opRO = 2*Nfft; 
    N_opR = N_opR + Nfft*2; % bit flipping
end

N_opM*n_samp/(Nfft-ParamRxDSP.hilbert_tap)
N_opA*n_samp/(Nfft-ParamRxDSP.hilbert_tap)


ParamECA_DSP.E_KK =  n_samp*(N_opM*ParamCMOS.E_opM+N_opA*ParamCMOS.E_opA+...
    N_opRO*n_ADC*ParamCMOS.E_opRO+N_opR*ParamCMOS.E_opR+N_opG*ParamCMOS.E_opG*n_ADC)...
    /(R_code*(Nfft-ParamRxDSP.hilbert_tap)*log2(MM));

%% Rx-DSP Resample at Rx after KK (downconversion + LPF+CDC)
N_h = ParamRxDSP.CD_tap;
n_parallel = 2*ParamSig.Baud_Rate/f_clk;
if ParamRxDSP.KKoverSamp == 64/28
    n_samp = 2;
    Nfft1 = 1024;

    [N_opM1,N_opA1] = numOP_SplitRadix_FFT(Nfft1);
    
    N_opM1 = N_opM1/2;
    N_opA1 = N_opA1/2 + Nfft1*4/2;
    

    Nfft2 = 7;
    Nfft3 = 128;

    N_opM2 = 8*2;
    N_opA2 = 36*2;

    [N_opM3,N_opA3] = numOP_SplitRadix_FFT(Nfft3);

    N_opM = N_opM1 + Nfft2*N_opM3 + Nfft3*N_opM2;
    N_opA = N_opA1 + Nfft2*N_opA3 + Nfft3*N_opA2;
    Nfft_out = Nfft2*Nfft3;
    Nfft_in = Nfft1;
elseif ParamRxDSP.KKoverSamp == 72/28
    %1152 -> 896
    n_samp = 2;
    Nfft1 = 9;
    Nfft2 = 128;
    
    N_opM1 = 20; 
    N_opA1 = 84;
    [N_opM2,N_opA2] = numOP_SplitRadix_FFT(Nfft2);
    N_opM = Nfft2*N_opM1 + Nfft1*N_opM2;
    N_opA = Nfft2*N_opA1 + Nfft1*N_opA2;
   
    N_opM = N_opM/2;
    N_opA = N_opA/2 + Nfft1*Nfft2*4/2;


    Nfft3 = 7;
    Nfft4 = 128;

    N_opM3 = 8*2;
    N_opA3 = 36*2;

    [N_opM4,N_opA4] = numOP_SplitRadix_FFT(Nfft4);

    N_opM = N_opM + Nfft3*N_opM4 + Nfft4*N_opM3;
    N_opA = N_opA + Nfft3*N_opA4 + Nfft4*N_opA3;
    
    
    
    Nfft_out = Nfft3*Nfft4;
    Nfft_in = Nfft1*Nfft2;

elseif ParamRxDSP.KKoverSamp == 80/28
    %1280 -> 896
    n_samp = 2;
    Nfft1 = 5;
    Nfft2 = 256;
    
    N_opM1 = 5*2; 
    N_opA1 = 17*2;
    [N_opM2,N_opA2] = numOP_SplitRadix_FFT(Nfft2);
    N_opM = Nfft2*N_opM1 + Nfft1*N_opM2;
    N_opA = Nfft2*N_opA1 + Nfft1*N_opA2;
   
    N_opM = N_opM/2;
    N_opA = N_opA/2 + Nfft1*Nfft2*4/2;


    Nfft3 = 7;
    Nfft4 = 128;

    N_opM3 = 8*2;
    N_opA3 = 36*2;

    [N_opM4,N_opA4] = numOP_SplitRadix_FFT(Nfft4);

    N_opM = N_opM + Nfft3*N_opM4 + Nfft4*N_opM3;
    N_opA = N_opA + Nfft3*N_opA4 + Nfft4*N_opA3;
    
    
    
    Nfft_out = Nfft3*Nfft4;
    Nfft_in = Nfft1*Nfft2;
elseif ParamRxDSP.KKoverSamp == 88/28
    %1408 -> 896
    n_samp = 2;
    Nfft1 = 11;
    Nfft2 = 128;
    
    N_opM1 = 40; 
    N_opA1 = 168;
    [N_opM2,N_opA2] = numOP_SplitRadix_FFT(Nfft2);
    N_opM = Nfft2*N_opM1 + Nfft1*N_opM2;
    N_opA = Nfft2*N_opA1 + Nfft1*N_opA2;
   
    N_opM = N_opM/2;
    N_opA = N_opA/2 + Nfft1*Nfft2*4/2;


    Nfft3 = 7;
    Nfft4 = 128;

    N_opM3 = 8*2;
    N_opA3 = 36*2;

    [N_opM4,N_opA4] = numOP_SplitRadix_FFT(Nfft4);

    N_opM = N_opM + Nfft3*N_opM4 + Nfft4*N_opM3;
    N_opA = N_opA + Nfft3*N_opA4 + Nfft4*N_opA3;
    
    
    
    Nfft_out = Nfft3*Nfft4;
    Nfft_in = Nfft1*Nfft2;
end

% pointwise multiplication

N_tap = ceil(Nfft_out/2*(1+ParamSig.roll_off));

N_opM = N_opM + 3*N_tap;
N_opA = N_opA + 3*N_tap;


N_opM*n_samp/(Nfft_out-ParamRxDSP.CD_tap)
N_opA*n_samp/(Nfft_out-ParamRxDSP.CD_tap)

N_opR = 5*Nfft_in+4*(Nfft_out-N_h);

N_opG = (Nfft_in+4*(Nfft_out-N_h))*ceil(log2(n_parallel));

ParamECA_DSP.E_Rx_resamp = n_samp*(N_opM*ParamCMOS.E_opM+N_opA*ParamCMOS.E_opA+n_ADC*N_opR*ParamCMOS.E_opR+N_opG*ParamCMOS.E_opG*n_ADC)/...
    ((Nfft_out-ParamRxDSP.CD_tap)*R_code*log2(MM));



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
N_opM = 3*(2*n_tap);
N_opA = 5*(2*n_tap)  + 2*(n_tap-1) + 2 + 2*n_tap;
N_opR = 2*2*(n_tap + n_tap );
N_opG = 0;
ParamECA_DSP.E_FFE_adpt = n_samp*(N_opM*ParamCMOS.E_opM+N_opA*ParamCMOS.E_opA+N_opR*ParamCMOS.E_opR...
       +N_opG*ParamCMOS.E_opG)/(R_code*log2(MM));

%% Convert DSP Ener-consump to Pwr-consump
ParamSig.Info_Bit_Rate = 100e9; 
ParamPCA.P_RC = ParamECA_DSP.E_RC*1e-15*ParamSig.Info_Bit_Rate;
ParamPCA.P_Tx_resamp = ParamECA_DSP.E_Tx_resamp*1e-15*ParamSig.Info_Bit_Rate; % complex
ParamPCA.P_KK = ParamECA_DSP.E_KK*1e-15*ParamSig.Info_Bit_Rate;
ParamPCA.P_Rx_resamp = ParamECA_DSP.E_Rx_resamp*1e-15*ParamSig.Info_Bit_Rate;
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


ParamPCA.P_DAC = num_DAC*560e-3/60e9*ParamDAC.DAC_Rate/8*ParamDAC.qnbit_DAC ;
ParamPCA.P_ADC = num_ADC*303e-3/60e9*ParamSig.Baud_Rate*ParamRxDSP.KKoverSamp/7*ParamADC.qnbit_ADC;

if ParamOpt.BER_target == 3.8e-3
    ParamPCA.P_FEC = (0.42*0.05+0.42*0.95*(ParamSig.Info_Bit_Rate/400e9))/2;
else
    ParamPCA.P_FEC = 0.42/4;
end
ParamPCA.P_PD = 0.169 ;
...+ 10^(ParamPhysicalModel.Received_Pwr_dBm/10)/1e3*0.76*5;

% ParamPhysicalModel.ModVpp = ParamSys.Vpi*ParamSys.Vpp_over_Vpi;
% if ParamControl.VSB_or_Not
%     ParamPhysicalModel.ModVrms = ParamPhysicalModel.Vrms_over_VM*ParamPhysicalModel.ModVpp/2;
% else
%     ParamPhysicalModel.ModVrms = mean([ParamPhysicalModel.Vrms_over_VM_I,ParamPhysicalModel.Vrms_over_VM_Q])*ParamPhysicalModel.ModVpp/2;
% end
% 
% ParamPhysicalModel.ModVbias = ParamSys.Vpi*ParamSys.Vbias_over_Vpi;
% % ParamPCA.P_mod = num_Mod*((ParamPhysicalModel.ModVrms)^2/50 + (ParamPhysicalModel.ModVbias)^2/50 );
% ParamPCA.P_mod = num_Mod*(ParamPhysicalModel.ModVrms)^2/50;


ParamPCA.P_driver = num_driver*2.7/4;

ParamPCA.P_laser_TEC = 1.2*2.4;
switch ParamControl.Laser_case
    case 1
        ParamPCA.P_laser = 0.08*1.8;
    case 2
        ParamPCA.P_laser = 0.15*1.8;
    case 3
        ParamPCA.P_laser = ParamLas.current/1000*ParamLas.Voltage;
end


P_all = [ParamPCA.P_DAC,ParamPCA.P_ADC,...
         ParamPCA.P_RC,ParamPCA.P_Tx_resamp,...  
         ParamPCA.P_KK,...
         ParamPCA.P_Rx_resamp, ...
         ParamPCA.P_TR, ParamPCA.P_FFE_adpt,...
         ParamPCA.P_FEC,...
         ParamPCA.P_laser_TEC,ParamPCA.P_laser...
         ParamPCA.P_driver,...
         ParamPCA.P_PD];
ParamPCA.P_sum = sum(P_all);
% disp(ParamPCA);

if ParamControl.Plot_Power_Pie_or_Not
    label_all = {'DAC','ADC','RC','Tx-Freq','KK',...
                 'Rx-Freq','TR','FFE-adapt',...
                 'FEC','laser(TEC)','laser(tunable)','modulator+driver','PIN-TIA'};
    explode = ones(size(P_all));
    figure; pie(P_all,explode,label_all);
    

end



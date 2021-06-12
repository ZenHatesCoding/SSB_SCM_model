PCA_Param_Init;

FFT_sharing_or_Not = 1;
%% KK
n_samp = 2;
N_FIR = 28+1;
% hilbert transform
Nfft = 1024;

[N_opM,N_opA] = numOP_SplitRadix_FFT(Nfft);
if FFT_sharing_or_Not == 0
    
    % single-tap mul + FFT/IFFT
    % multiplied by ii/-1i, swtich real and imag
    N_opM_iterKK = 3*N_opM + 2*3*Nfft;
    N_opM_cKK = 2*N_opM;

    N_opA_iterKK = 3*N_opA + 2*3*Nfft;
    N_opA_cKK = 2*N_opA;
else
    
    N_opM_iterKK = N_opM/2 + N_opM + 2*3*Nfft;
    N_opM_cKK = N_opM;
    
    N_opA_iterKK = N_opA/2 + N_opA + 2*3*Nfft + 4*Nfft;
    N_opA_cKK = N_opA;
    
end

N_opM_iterKK = N_opM_iterKK*n_samp/(Nfft-N_FIR+1);
N_opM_cKK  = N_opM_cKK*n_samp/Nfft;

N_opA_iterKK = N_opA_iterKK*n_samp/(Nfft-N_FIR+1);
N_opA_cKK  = N_opA_cKK*n_samp/Nfft;

E_iterKK = N_opM_iterKK*ParamCMOS.E_opM+N_opA_iterKK*ParamCMOS.E_opA;
E_cKK = N_opM_cKK*ParamCMOS.E_opM+N_opA_cKK*ParamCMOS.E_opA;
%% Volterra FFE
L2 = 3;
M2 = 2*L2+1; 
N_opM_Vol = 2*M2*(M2+1)/2;
N_opA_Vol = M2*(M2+1)/2;

N_opM_Vol_adpt = 3*M2*(M2+1)/2;
N_opA_Vol_adpt = 2*M2*(M2+1)/2+1;

E_cKK_fix = E_cKK + N_opM_Vol*ParamCMOS.E_opM+N_opA_Vol*ParamCMOS.E_opA;
E_cKK_adpt = E_cKK + N_opM_Vol_adpt*ParamCMOS.E_opM+N_opA_Vol_adpt*ParamCMOS.E_opA;
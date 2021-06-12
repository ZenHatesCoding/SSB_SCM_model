function [N_opM, N_opA] = numOP_SplitRadix_FFT(Nfft)
    N_opM = Nfft*log2(Nfft)-3*Nfft+4;
    N_opA = 3*Nfft*log2(Nfft)-3*Nfft+4;
end
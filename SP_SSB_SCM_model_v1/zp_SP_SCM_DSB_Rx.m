function [BER_avg, BER_list, SNR_list] = zp_SP_SCM_DSB_Rx(ParamControl,ParamRxDSP,ParamSig,ParamADC,ParamPS,ParamFib,ParamMod,ParamPD,Rx_I,Tx_Sequence)    
    
    BER_list = zeros(size(ParamSig.SC));
    SNR_list = zeros(size(ParamSig.SC));
    
    
    %% Plot spectrum after PD
    for idx = 1:length(ParamSig.SC)
        Tx_Sequence_SC = Tx_Sequence{idx};
   
        [BER,SNR] = zp_SP_SC_PAMx_DSB_Rx(ParamControl,ParamRxDSP,ParamSig,ParamADC,ParamPS,ParamFib,Rx_I,Tx_Sequence_SC,idx,ParamSig.SC(idx));

        BER_list(idx) = BER;
        SNR_list(idx) = SNR;
    end   
    BER_avg = BER_list.*ParamSig.SC_Baud_Rate.*ParamSig.SC./sum(ParamSig.SC_Baud_Rate.*ParamSig.SC);
end
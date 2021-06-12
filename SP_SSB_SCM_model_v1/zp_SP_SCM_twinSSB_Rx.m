function [BER_avg, BER_list, SNR_list] = zp_SP_SCM_twinSSB_Rx(ParamControl,ParamRxDSP,ParamSig,ParamADC,ParamPS1,ParamPS2,ParamFib,ParamMod,ParamPD,ParamVSB,Rx_I1,Rx_I2,Tx_Sequence1,Tx_Sequence2)    
    
    BER_list = zeros(size(ParamSig.SC));
    SNR_list = zeros(size(ParamSig.SC));
    
    
    %% Plot spectrum after PD
    for idx = 1:length(ParamSig.SC)
        Tx_Sequence_SC1 = Tx_Sequence1{idx};
        Tx_Sequence_SC2 = Tx_Sequence2{idx};
        [BER,SNR] = zp_SP_SC_QAMx_twinSSB_Rx(ParamControl,ParamRxDSP,ParamSig,ParamADC,ParamPS1,ParamPS2,ParamFib,ParamMod,ParamPD,ParamVSB,Rx_I1,Rx_I2,Tx_Sequence_SC1,Tx_Sequence_SC2,idx,ParamSig.SC(idx));

        BER_list(idx) = BER;
        SNR_list(idx) = SNR;
    end   
    BER_avg = BER_list.*ParamSig.SC_Baud_Rate.*ParamSig.SC./sum(ParamSig.SC_Baud_Rate.*ParamSig.SC);
end
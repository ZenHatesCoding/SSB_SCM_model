%% Update by ZP on 02/06/2019
% OFC 2018 W4E.2
% ***************************************************
% achieved using hilbert transform
% ***************************************************
function Rx_Time_Data = rx_KK_woUpSamp_detection(Rx_I,ParamControl)
    sqrt_Rx_I = sqrt(Rx_I);
    E0 = sqrt(mean(Rx_I));
    Rx_Time_Data = sqrt_Rx_I-...
                   E0/2*imag(hilbert(real(2*sqrt_Rx_I/E0-1/2*Rx_I/(E0^2)))).^2;
    if ParamControl.KK_real_output_or_Not == 0
        Rx_Time_Data = Rx_Time_Data + 1i*imag(hilbert(Rx_Time_Data));
    end
end
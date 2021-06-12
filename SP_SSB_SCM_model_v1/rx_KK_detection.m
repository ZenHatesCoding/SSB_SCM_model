% ***************************************************
% achieved using hilbert transform
% ***************************************************
function Rx_Time_Data = rx_KK_detection(Rx_I,ParamControl)
    % trick for AC coupling PD
%     alpha = 1.7;
%     Rx_I = Rx_I - mean(Rx_I);
%     Rx_I = Rx_I + alpha*abs(min(Rx_I));
    
    log_Rx_I = log(abs(Rx_I));
    Rx_Phi =  1/2 * imag(hilbert(real(log_Rx_I)));
    if ParamControl.KK_real_output_or_Not == 1
        Rx_Time_Data = sqrt(Rx_I).*cos(Rx_Phi);
    else
        Rx_Time_Data = sqrt(Rx_I).*exp(1j*Rx_Phi);
    end
    
end
% ***************************************************
% achieved using 32-tap FIR filter
% ***************************************************
function Rx_Time_Data = rx_KK_detection_FIR_hilbert(Rx_I,N)

    sqrt_Rx_I = sqrt(abs(Rx_I));
    log_sqrt_Rx_I = log(sqrt_Rx_I);
    h = zeros(1,N);
    for n = -(N-1)/2:(N-1)/2
        if n ~=0
            h(n+(N-1)/2+1) = 2/pi*(sin(pi*n/2))^2/n;
        end
    end
    hilbert_log_sqrt_Rx_I = conv(log_sqrt_Rx_I,h,'same');    
    Rx_Time_Data = sqrt_Rx_I.*exp(1i*hilbert_log_sqrt_Rx_I);
    
end
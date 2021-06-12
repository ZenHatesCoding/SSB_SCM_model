% v2: to be used after cKK
% ***************************************************
% achieved using hilbert transform
% ***************************************************
function Rx_Time_Data = rx_iter_KK_detection_v2(Rx_I,E_r0,iter)

    
    E0 = mean(sqrt(Rx_I));
    E_r = real(E_r0 - mean(E_r0));
    
    for i = 1:iter        
        E = hilbert(real(E_r));
        Ef = fftshift(fft(E));
        Ef(1:floor(length(Ef)/2)) = 0;
        E = ifft(ifftshift(Ef));
        E_i = imag(E);
        E_r = sqrt(abs(Rx_I-E_i.^2))-E0;
    end
%     Rx_Time_Data = E_r + 1i*E_i;
    Rx_Time_Data = E_r;
    
end
function out_signal = Add_AWGN(signal,SNR)

    signal_power = (norm(signal,2))^2/(length(signal)); 
    % norm(x,2) = sqrt(sum(x.^2));
%     noise_power = 10*log10((signal_power/(10^(SNR/10)))*IFFT_bin_length/carrier_count);
    noise_power = 10*log10((signal_power/(10^(SNR/10))));
    
    noise = wgn(1,length(signal),noise_power,'complex');
    out_signal = signal+noise;

end


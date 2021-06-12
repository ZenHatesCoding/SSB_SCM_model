function dec_symbols = rx_16QAM_Decision(rx_symbols)



n = length(rx_symbols);

dec_symbols = zeros(1,n); 

% boundaries: 3/sqrt(10) = 0.9487 2/sqrt(10) = 0.6325 1/sqrt(10) = 0.3162 
for i = 1:n
    if real(rx_symbols(i)) >= 0.6325
        dec_symbols(i) = dec_symbols(i)+0.9487;
    elseif real(rx_symbols(i)) >= 0 && real(rx_symbols(i)) < 0.6325
        dec_symbols(i) = dec_symbols(i)+0.3162;
    elseif real(rx_symbols(i)) >= -0.6325 && real(rx_symbols(i)) < 0
        dec_symbols(i) = dec_symbols(i)-0.3162;
    else 
        dec_symbols(i) = dec_symbols(i)-0.9487; 
    end
    
    if imag(rx_symbols(i)) >= 0.6325
        dec_symbols(i) = dec_symbols(i)+1j*0.9487;
    elseif imag(rx_symbols(i)) >= 0 && imag(rx_symbols(i)) < 0.6325
        dec_symbols(i) = dec_symbols(i)+1j*0.3162;
    elseif imag(rx_symbols(i)) >= -0.6325 && imag(rx_symbols(i)) < 0
        dec_symbols(i) = dec_symbols(i)-1j*0.3162;
    else 
        dec_symbols(i) = dec_symbols(i)-1j*0.9487;
    end
end
end


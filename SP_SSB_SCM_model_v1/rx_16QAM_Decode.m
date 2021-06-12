function rx_bits = rx_16QAM_Decode(rx_symbols)



n = length(rx_symbols);

rx_bits = zeros(1,4*n);

for i = 1:n
    if real(rx_symbols(i)) >= 0.6325
        rx_bits(4*(i-1)+1) = 1;
        rx_bits(4*(i-1)+2) = 0;
    elseif real(rx_symbols(i)) >= 0 && real(rx_symbols(i)) < 0.6325
        rx_bits(4*(i-1)+1) = 1;
        rx_bits(4*(i-1)+2) = 1;
    elseif real(rx_symbols(i)) >= -0.6325 && real(rx_symbols(i)) < 0
        rx_bits(4*(i-1)+1) = 0;
        rx_bits(4*(i-1)+2) = 1;
    else 
        rx_bits(4*(i-1)+1) = 0;
        rx_bits(4*(i-1)+2) = 0;
    end
    
    if imag(rx_symbols(i)) >= 0.6325
        rx_bits(4*(i-1)+3) = 1;
        rx_bits(4*(i-1)+4) = 0;
    elseif imag(rx_symbols(i)) >= 0 && imag(rx_symbols(i)) < 0.6325
        rx_bits(4*(i-1)+3) = 1;
        rx_bits(4*(i-1)+4) = 1;
    elseif imag(rx_symbols(i)) >= -0.6325 && imag(rx_symbols(i)) < 0
        rx_bits(4*(i-1)+3) = 0;
        rx_bits(4*(i-1)+4) = 1;
    else 
        rx_bits(4*(i-1)+3) = 0;
        rx_bits(4*(i-1)+4) = 0;
    end
end
end


function dec_symbols = rx_64QAM_Decision(rx_symbols)
    n = length(rx_symbols);
    dec_symbols = zeros(1,n); 
    % sum(norm(table)^2)/length(table) = 42
    % boundaries: x/sqrt(42) 
    for i = 1:n
        if real(rx_symbols(i)) >= 6/sqrt(42)
            dec_symbols(i) = 7/sqrt(42);
        elseif real(rx_symbols(i)) >= 4/sqrt(42) && real(rx_symbols(i)) < 6/sqrt(42)
            dec_symbols(i) = 5/sqrt(42);
        elseif real(rx_symbols(i)) >= 2/sqrt(42) && real(rx_symbols(i)) < 4/sqrt(42)
            dec_symbols(i) = 3/sqrt(42);
        elseif real(rx_symbols(i)) >= 0 && real(rx_symbols(i)) < 2/sqrt(42)
            dec_symbols(i) = 1/sqrt(42);
        elseif real(rx_symbols(i)) >= -2/sqrt(42) && real(rx_symbols(i)) < 0
            dec_symbols(i) = -1/sqrt(42);
        elseif real(rx_symbols(i)) >= -4/sqrt(42) && real(rx_symbols(i)) < -2/sqrt(42)
            dec_symbols(i) = -3/sqrt(42);   
        elseif real(rx_symbols(i)) >= -6/sqrt(42) && real(rx_symbols(i)) < -4/sqrt(42)
            dec_symbols(i) = -5/sqrt(42);   
        else 
            dec_symbols(i) = -7/sqrt(42); 
        end

        if imag(rx_symbols(i)) >= 6/sqrt(42)
            dec_symbols(i) = dec_symbols(i)+1j*7/sqrt(42);
        elseif imag(rx_symbols(i)) >= 4/sqrt(42) && imag(rx_symbols(i)) < 6/sqrt(42)
            dec_symbols(i) = dec_symbols(i)+1j*5/sqrt(42);
        elseif imag(rx_symbols(i)) >= 2/sqrt(42) && imag(rx_symbols(i)) < 4/sqrt(42)
            dec_symbols(i) = dec_symbols(i)+1j*3/sqrt(42);
        elseif imag(rx_symbols(i)) >= 0 && imag(rx_symbols(i)) < 2/sqrt(42)
            dec_symbols(i) = dec_symbols(i)+1j*1/sqrt(42);
        elseif imag(rx_symbols(i)) >= -2/sqrt(42) && imag(rx_symbols(i)) < 0
            dec_symbols(i) = dec_symbols(i)+1j*(-1/sqrt(42));
        elseif imag(rx_symbols(i)) >= -4/sqrt(42) && imag(rx_symbols(i)) < -2/sqrt(42)
            dec_symbols(i) = dec_symbols(i)+1j*(-3/sqrt(42));   
        elseif imag(rx_symbols(i)) >= -6/sqrt(42) && imag(rx_symbols(i)) < -4/sqrt(42)
            dec_symbols(i) = dec_symbols(i)+1j*(-5/sqrt(42));   
        else 
            dec_symbols(i) = dec_symbols(i)+1j*(-7/sqrt(42)); 
        end
        
    end
end
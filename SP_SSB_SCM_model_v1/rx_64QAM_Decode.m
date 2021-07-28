function rx_bits = rx_64QAM_Decode(rx_symbols)
    n = length(rx_symbols);
    rx_bits = zeros(1,6*n);

    % sum(norm(table)^2)/length(table) = 42
    % boundaries: x/sqrt(42) 
    
    for i = 1:n
        if real(rx_symbols(i)) >= 6/sqrt(42)
            rx_bits(6*(i-1)+(1:3)) = [1 0 0];
        elseif real(rx_symbols(i)) >= 4/sqrt(42) && real(rx_symbols(i)) < 6/sqrt(42)
            rx_bits(6*(i-1)+(1:3)) = [1 0 1];
        elseif real(rx_symbols(i)) >= 2/sqrt(42) && real(rx_symbols(i)) < 4/sqrt(42)
            rx_bits(6*(i-1)+(1:3)) = [1 1 1];
        elseif real(rx_symbols(i)) >= 0 && real(rx_symbols(i)) < 2/sqrt(42)
            rx_bits(6*(i-1)+(1:3)) = [1 1 0];
        elseif real(rx_symbols(i)) >= -2/sqrt(42) && real(rx_symbols(i)) < 0
            rx_bits(6*(i-1)+(1:3)) = [0 1 0];
        elseif real(rx_symbols(i)) >= -4/sqrt(42) && real(rx_symbols(i)) < -2/sqrt(42)
            rx_bits(6*(i-1)+(1:3)) = [0 1 1];  
        elseif real(rx_symbols(i)) >= -6/sqrt(42) && real(rx_symbols(i)) < -4/sqrt(42)
            rx_bits(6*(i-1)+(1:3)) = [0 0 1];
        else 
            rx_bits(6*(i-1)+(1:3)) = [0 0 0];
        end
        
       
        
        if imag(rx_symbols(i)) >= 6/sqrt(42)
            rx_bits(6*(i-1)+(4:6)) = [1 0 0];
        elseif imag(rx_symbols(i)) >= 4/sqrt(42) && imag(rx_symbols(i)) < 6/sqrt(42)
            rx_bits(6*(i-1)+(4:6)) = [1 0 1];
        elseif imag(rx_symbols(i)) >= 2/sqrt(42) && imag(rx_symbols(i)) < 4/sqrt(42)
            rx_bits(6*(i-1)+(4:6)) = [1 1 1];
        elseif imag(rx_symbols(i)) >= 0 && imag(rx_symbols(i)) < 2/sqrt(42)
            rx_bits(6*(i-1)+(4:6)) = [1 1 0];
        elseif imag(rx_symbols(i)) >= -2/sqrt(42) && imag(rx_symbols(i)) < 0
            rx_bits(6*(i-1)+(4:6)) = [0 1 0];
        elseif imag(rx_symbols(i)) >= -4/sqrt(42) && imag(rx_symbols(i)) < -2/sqrt(42)
            rx_bits(6*(i-1)+(4:6)) = [0 1 1];  
        elseif imag(rx_symbols(i)) >= -6/sqrt(42) && imag(rx_symbols(i)) < -4/sqrt(42)
            rx_bits(6*(i-1)+(4:6)) = [0 0 1];   
        else 
            rx_bits(6*(i-1)+(4:6)) = [0 0 0];
        end
        
    end
end


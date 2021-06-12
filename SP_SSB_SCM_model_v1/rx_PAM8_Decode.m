function rx_bits = rx_PAM8_Decode(rx_symbols)



n = length(rx_symbols);

rx_bits = zeros(1,3*n);

% PAM8 [0  1  3  2  6  7  5  4]
%      -7 -5 -3 -1  1  3  5  7 
table = -7:2:7;
norm_coef = sqrt(norm(table)^2/length(table));

for i = 1:n
    if real(rx_symbols(i)) >= 6/norm_coef
        rx_bits(3*(i-1)+1) = 1;
        rx_bits(3*(i-1)+2) = 0;
        rx_bits(3*(i-1)+3) = 0;
    elseif real(rx_symbols(i)) >= 4/norm_coef && real(rx_symbols(i)) < 6/norm_coef
        rx_bits(3*(i-1)+1) = 1;
        rx_bits(3*(i-1)+2) = 0;
        rx_bits(3*(i-1)+3) = 1;
    elseif real(rx_symbols(i)) >= 2/norm_coef && real(rx_symbols(i)) < 4/norm_coef
        rx_bits(3*(i-1)+1) = 1;
        rx_bits(3*(i-1)+2) = 1;
        rx_bits(3*(i-1)+3) = 1;
    elseif real(rx_symbols(i)) >= 0/norm_coef && real(rx_symbols(i)) < 2/norm_coef
        rx_bits(3*(i-1)+1) = 1;
        rx_bits(3*(i-1)+2) = 1;
        rx_bits(3*(i-1)+3) = 0;
    elseif real(rx_symbols(i)) >= -2/norm_coef && real(rx_symbols(i)) < 0/norm_coef
        rx_bits(3*(i-1)+1) = 0;
        rx_bits(3*(i-1)+2) = 1;
        rx_bits(3*(i-1)+3) = 0;
    elseif real(rx_symbols(i)) >= -4/norm_coef && real(rx_symbols(i)) < -2/norm_coef
        rx_bits(3*(i-1)+1) = 0;
        rx_bits(3*(i-1)+2) = 1;
        rx_bits(3*(i-1)+3) = 1;
    elseif real(rx_symbols(i)) >= -6/norm_coef && real(rx_symbols(i)) < -4/norm_coef
        rx_bits(3*(i-1)+1) = 0;
        rx_bits(3*(i-1)+2) = 0;
        rx_bits(3*(i-1)+3) = 1;        
    else 
        rx_bits(3*(i-1)+1) = 0;
        rx_bits(3*(i-1)+2) = 0;
        rx_bits(3*(i-1)+3) = 0; 
    end

end
end


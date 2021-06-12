
function rx_bits = rx_8QAM_Decode(rx_symbols)

    n = length(rx_symbols);

    rx_bits = zeros(1,3*n);

    for i = 1:n

        table = [-1+1j, 1+1j, -1-1j, 1-1j, -1-sqrt(3), 1j*(1+sqrt(3)), -1j*(1+sqrt(3)), 1+sqrt(3)]/sqrt((3+sqrt(3)));
    %     table = [-1-1j,  -1j*(1+sqrt(3)), 1+sqrt(3), 1-1j, -1-sqrt(3), -1+1j, 1+1j, 1j*(1+sqrt(3))]/sqrt((3+sqrt(3)));
%         table = pwr_normalization(table);
        A = abs(rx_symbols(i) - table);
        [B I] = min(A);
        rx_bits(3*(i-1)+1) = floor((I-1)/4);
        rx_bits(3*(i-1)+2) = floor((I-rx_bits(3*(i-1)+1)*4-1)/2);
        rx_bits(3*(i-1)+3) = floor((I-rx_bits(3*(i-1)+1)*4-rx_bits(3*(i-1)+2)*2-1)/1);
    end
end
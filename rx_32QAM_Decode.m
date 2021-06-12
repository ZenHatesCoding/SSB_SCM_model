function rx_bits = rx_32QAM_Decode(rx_symbols)



n = length(rx_symbols);

rx_bits = zeros(1,5*n);

table = [-3-5*1j -1-5*1j -3+5*1j -1+5*1j -5-3*1j -5-1*1j -5+3*1j -5+1*1j...
        -1-3*1j -1-1*1j -1+3*1j -1+1*1j -3-3*1j -3-1*1j -3+3*1j -3+1*1j...
        3-5*1j 1-5*1j 3+5*1j 1+5*1j 5-3*1j 5-1*1j 5+3*1j 5+1*1j...
        1-3*1j 1-1*1j 1+3*1j 1+1*1j 3-3*1j 3-1*1j 3+3*1j 3+1*1j];
%      table = [-7-3*1j -7-1*1j -7+3*1j -7+1*1j -5-3*1j -5-1*1j -5+3*1j -5+1*1j...
%         -1-3*1j -1-1*1j -1+3*1j -1+1*1j -3-3*1j -3-1*1j -3+3*1j -3+1*1j...
%         7-3*1j 7-1*1j 7+3*1j 7+1*1j 5-3*1j 5-1*1j 5+3*1j 5+1*1j...
%         1-3*1j 1-1*1j 1+3*1j 1+1*1j 3-3*1j 3-1*1j 3+3*1j 3+1*1j]; 
table = table/sqrt(mean(abs(table).^2));

for i = 1:n
    
    A = abs(rx_symbols(i) - table);
    [B I] = min(A);
    
%     AA = dec2bin(I,5);
%     for bb = 1:5
%         rx_bits(5*(i-1)+bb) = str2num(AA(bb));
% 
%     end
    
    rx_bits(5*(i-1)+1) = floor((I-1)/16);
    rx_bits(5*(i-1)+2) = floor((I-rx_bits(5*(i-1)+1)*16-1)/8);
    rx_bits(5*(i-1)+3) = floor((I-rx_bits(5*(i-1)+1)*16-rx_bits(5*(i-1)+2)*8-1)/4);
    rx_bits(5*(i-1)+4) = floor((I-rx_bits(5*(i-1)+1)*16-rx_bits(5*(i-1)+2)*8-rx_bits(5*(i-1)+3)*4-1)/2);
    rx_bits(5*(i-1)+5) = floor((I-rx_bits(5*(i-1)+1)*16-rx_bits(5*(i-1)+2)*8-rx_bits(5*(i-1)+3)*4-rx_bits(5*(i-1)+4)*2-1)/1);
end
end


function dec_symbols = rx_32QAM_Decision(rx_symbols)



n = length(rx_symbols);

dec_symbols = zeros(1,n);

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
    
    dec_symbols(i) = table(I); 
end
end


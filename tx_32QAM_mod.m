function  Mod_Symbols = tx_32QAM_mod(Bits_In)

    L = length(Bits_In);

%     00000 00001 00010 00011 00100 00101 00110 00111
%     01000 01001 01010 01011 01100 01101 01110 01111
%     10000 10001 10010 10011 10100 10101 10110 10111
%     11000 11001 11010 11011 11100 11101 11110 11111
    table = [-3-5*1j -1-5*1j -3+5*1j -1+5*1j -5-3*1j -5-1*1j -5+3*1j -5+1*1j...
        -1-3*1j -1-1*1j -1+3*1j -1+1*1j -3-3*1j -3-1*1j -3+3*1j -3+1*1j...
        3-5*1j 1-5*1j 3+5*1j 1+5*1j 5-3*1j 5-1*1j 5+3*1j 5+1*1j...
        1-3*1j 1-1*1j 1+3*1j 1+1*1j 3-3*1j 3-1*1j 3+3*1j 3+1*1j];

%      table = [-7-3*1j -7-1*1j -7+3*1j -7+1*1j -5-3*1j -5-1*1j -5+3*1j -5+1*1j...
%         -1-3*1j -1-1*1j -1+3*1j -1+1*1j -3-3*1j -3-1*1j -3+3*1j -3+1*1j...
%         7-3*1j 7-1*1j 7+3*1j 7+1*1j 5-3*1j 5-1*1j 5+3*1j 5+1*1j...
%         1-3*1j 1-1*1j 1+3*1j 1+1*1j 3-3*1j 3-1*1j 3+3*1j 3+1*1j]; 
    
    
    table = table/sqrt(mean(abs(table).^2));
        

    inp=reshape(Bits_In,5,L/5);
    Mod_Symbols=table([16 8 4 2 1]*inp+1); 

end
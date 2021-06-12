function  Mod_Symbols = tx_8QAM_mod(Bits_In)

    L = length(Bits_In);

    table = [-1+1j, 1+1j, -1-1j, 1-1j, -1-sqrt(3), 1j*(1+sqrt(3)), -1j*(1+sqrt(3)), 1+sqrt(3)]/sqrt((3+sqrt(3)));
%     table = [-1-1j,  -1j*(1+sqrt(3)), 1+sqrt(3), 1-1j, -1-sqrt(3), -1+1j, 1+1j, 1j*(1+sqrt(3))]/sqrt((3+sqrt(3)));
    
    table = pwr_normalization(table);
    inp=reshape(Bits_In,3,L/3);
    Mod_Symbols=table([4 2 1]*inp+1); 
    
end
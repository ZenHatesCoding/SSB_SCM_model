function  mod_symbols = tx_PAM8_mod(bits_in)
   full_len = length(bits_in);
   table = -7:2:7;
   norm_coef = sqrt(norm(table)^2/length(table));
   table = table/norm_coef; 
   table = table([0 1 3 2 6 7 5 4]+1); % Gray Code mapping pattern

   inp=reshape(bits_in,3,full_len/3); % PAM4 3 bits/symbol
                                      % 3 row full_len/3 col, #col =
                                      % #symbol
   
   mod_symbols=table([4 2 1]*inp+1);  % maps transmitted bits into PAM4 
end 
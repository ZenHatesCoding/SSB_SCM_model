function  mod_symbols = tx_16QAM_mod(bits_in)
full_len = length(bits_in);
   m=1;
   for k=-3:2:3
      for l=-3:2:3
         table(m) = (k+j*l)/sqrt(10); % power normalization
         m=m+1;
      end;
   end;
   table=table([0 1 3 2 4 5 7 6 12 13 15 14 8 9 11 10]+1); % Gray code mapping pattern for 8-PSK symbols
   inp=reshape(bits_in,4,full_len/4);
   mod_symbols=table([8 4 2 1]*inp+1);  % maps transmitted bits into 16QAM symbols
end 
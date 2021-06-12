function mod_data = tx_QPSK_mod(inf_bits)

   Bit_Length = length(inf_bits);

   table=exp(j*[1/4*pi -1/4*pi  3/4*pi -3/4*pi]);  % generates QPSK symbols
   inp=reshape(inf_bits,2,Bit_Length/2);
   mod_data=table([2 1]*inp+1);  % maps transmitted bits into QPSK symbols
end
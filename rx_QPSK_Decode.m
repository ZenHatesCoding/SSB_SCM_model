function rx_bits = rx_QPSK_Decode(rx_data)
   bit_len = length(rx_data)*2;
   re=real(rx_data);
   im=imag(rx_data);
   
   rx_bits(1:2:bit_len) = re<0;
   rx_bits(2:2:bit_len) = im<0;
end
function [rx_bits, num_parity_check_error] = rx_SP128_Decode(rx_symbols)
    n = length(rx_symbols);
    rx_bits = zeros(1,7*n/2);
    idx_2d_symbols = 1;
    
    dec_bound = [0.6325, 0, -0.6325,1i*0.6325,0,-1i*0.6325];
    num_parity_check_error = 0;
    for i = 1:2:n  
        rx_bits_1 = rx_16QAM_Decode(rx_symbols(i));
        rx_bits_2 = rx_16QAM_Decode(rx_symbols(i+1));
        
        
        if mod((sum(rx_bits_1)+sum(rx_bits_2)),2) == 0
            rx_bits(idx_2d_symbols:(idx_2d_symbols+3)) = rx_bits_1;
            rx_bits((idx_2d_symbols+4):(idx_2d_symbols+6)) = rx_bits_2(1:3);
        else
            num_parity_check_error = num_parity_check_error + 1;
            rx_Symbols_1_list = [real(rx_symbols(i)),real(rx_symbols(i)),real(rx_symbols(i)),...
                                 imag(rx_symbols(i)),imag(rx_symbols(i)),imag(rx_symbols(i))];
            rx_Symbols_2_list = [real(rx_symbols(i+1)),real(rx_symbols(i+1)),real(rx_symbols(i+1)),...
                                 imag(rx_symbols(i+1)),imag(rx_symbols(i+1)),imag(rx_symbols(i+1))];
            distance_1 = abs(rx_Symbols_1_list-dec_bound);
            distance_2 = abs(rx_Symbols_2_list-dec_bound);
            [min_distance_1,min_distance_idx_1] = min(distance_1);
            [min_distance_2,min_distance_idx_2] = min(distance_2);
            if min_distance_1 < min_distance_2
                if min_distance_idx_1 == 1
                    rx_bits_1(2) = flip_bit(rx_bits_1(2));
                elseif min_distance_idx_1 == 2
                    rx_bits_1(1) = flip_bit(rx_bits_1(1));
                elseif min_distance_idx_1 == 3
                    rx_bits_1(2) = flip_bit(rx_bits_1(2));
                elseif min_distance_idx_1 == 4
                    rx_bits_1(4) = flip_bit(rx_bits_1(4));
                elseif min_distance_idx_1 == 5
                    rx_bits_1(3) = flip_bit(rx_bits_1(3));
                elseif min_distance_idx_1 == 6
                    rx_bits_1(4) = flip_bit(rx_bits_1(4));
                end    
            else
                if min_distance_idx_2 == 1
                    rx_bits_2(2) = flip_bit(rx_bits_2(2));
                elseif min_distance_idx_2 == 2
                    rx_bits_2(1) = flip_bit(rx_bits_2(1));
                elseif min_distance_idx_2 == 3
                    rx_bits_2(2) = flip_bit(rx_bits_2(2));
                elseif min_distance_idx_2 == 4
                    rx_bits_2(4) = flip_bit(rx_bits_2(4));
                elseif min_distance_idx_2 == 5
                    rx_bits_2(3) = flip_bit(rx_bits_2(3));
                elseif min_distance_idx_2 == 6
                    rx_bits_2(4) = flip_bit(rx_bits_2(4));
                end  
            end
            rx_bits(idx_2d_symbols:(idx_2d_symbols+3)) = rx_bits_1;
            rx_bits((idx_2d_symbols+4):(idx_2d_symbols+6)) = rx_bits_2(1:3);
        end
        idx_2d_symbols = idx_2d_symbols+7;       
    end
end

function Xout = flip_bit(Xin)
    Xout = round(-(Xin-0.5)+0.5);
end


    
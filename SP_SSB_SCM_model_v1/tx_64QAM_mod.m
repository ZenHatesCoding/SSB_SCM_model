function mod_symbols = tx_64QAM_mod(bits_in)
    full_len = length(bits_in);
    % 4 12 28 20 52 60 44 36
    % 5 13 29 21 53 61 45 37
    % 7 15 31 23 55 63 47 39
    % 6 14 30 22 54 62 46 38
    % 2 10 26 18 50 58 42 34
    % 3 11 27 19 51 59 43 35
    % 1  9 25 17 49 57 41 33
    % 0  8 24 16 48 56 40 32
    
    m = 1;
    table = zeros(1,64);
    for k = -7:2:7
        for t = -7:2:7
            table(m) = 1i*k+t;
            m = m+1;
        end
    end
%     sum(norm(table)^2)/length(table) = 42
    table = pwr_normalization(table);
    gray_map = [ 0  8 24 16 48 56 40 32 ...
                 1  9 25 17 49 57 41 33 ...
                 3 11 27 19 51 59 43 35 ...
                 2 10 26 18 50 58 42 34 ...
                 6 14 30 22 54 62 46 38 ...
                 7 15 31 23 55 63 47 39 ...
                 5 13 29 21 53 61 45 37 ...
                 4 12 28 20 52 60 44 36];
    table = table(gray_map+1); % gray mapping
    inp=reshape(bits_in,6,full_len/6);
    mod_symbols=table([32 16 8 4 2 1]*inp+1);
    
    
end
% By ZP
% xn, received bit sequence,
% dn, desired bit sequence
% M number of taps
% mu update step
% Sps_in # Sample per symbol of input waveform 

function yn = QAMx_LMS_DD(xn, bn, SE, mu, Update_Tap_or_Not,Sps_in,LMS_Plot_Conv_or_Not)
    if size(xn,2)>1
        xn = xn.';
    end
    
    M = (length(bn)-1)/2;
    yn = zeros(round(length(xn)/Sps_in),1);
   
    y_idx = 1;
    
    if LMS_Plot_Conv_or_Not
        en_list = [];
        bn_list = [];
    end
    if Update_Tap_or_Not
        for n = (M+1):Sps_in:length(xn)-M
            y = bn.'*xn(n+M:-1:n-M);
%             y = bn.'*xn(n-M:1:n+M);
            yn(y_idx) = y;
            y_idx = y_idx + 1;
            if y_idx > length(yn)
                break;
            end
            switch SE
                case 5
                    en = rx_32QAM_Decision(y) - y;
                case 4
                    en = rx_16QAM_Decision(y) - y;
                case 3
                    en = tx_8QAM_mod(rx_8QAM_Decode(y)) - y;
                case 2
                    en = rx_QPSK_Decision(y) - y;
            end
            if LMS_Plot_Conv_or_Not
                en_list = [en_list,en];
                bn_list = [bn_list,bn(M)]; 
            end

%             bn = bn + mu*en*conj(xn(n-M:1:n+M));
            bn = bn + mu*en*conj(xn(n+M:-1:n-M));
        end
        if LMS_Plot_Conv_or_Not
            figure; plot(abs(en_list));
            figure; plot(abs(bn_list)); hold on; plot(angle(bn_list));
        end
    else
        for n = (M+1):Sps_in:length(xn)-M
            y = bn.'*xn(n+M:-1:n-M);
%             y = bn.'*xn(n-M:1:n+M);
            yn(y_idx) = y;
            y_idx = y_idx + 1;
            if y_idx > length(yn)
                break;
            end
        end
    end
    yn = yn.';
    
    
end
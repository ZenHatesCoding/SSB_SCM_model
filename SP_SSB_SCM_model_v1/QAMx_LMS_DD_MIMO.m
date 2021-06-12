%% Based on QAMx_LMS_DD
% Update by ZP 11/01/2020
% MIMO
% bnij : input j, output i
%%
% By ZP
% xn, received bit sequence,
% dn, desired bit sequence
% M number of taps
% mu update step
% Sps_in # Sample per symbol of input waveform 

function [yn1,yn2] = QAMx_LMS_DD_MIMO(xn1,xn2, bn11,bn12,bn21,bn22,...
    SE, mu, Update_Tap_or_Not,Sps_in,LMS_Plot_Conv_or_Not)
    if size(xn1,2)>1
        xn1 = xn1.';
    end
    if size(xn2,2)>1
        xn2 = xn2.';
    end
    
    M = (length(bn11)-1)/2;
    yn1 = zeros(round(length(xn1)/Sps_in),1);
    yn2 = zeros(round(length(xn2)/Sps_in),1);
    y_idx = 1;
    
    if LMS_Plot_Conv_or_Not
        en_list = [];
        bn_list = [];
    end
    
    for n = (M+1):Sps_in:length(xn1)-M
        y1 = bn11.'*xn1(n+M:-1:n-M)+bn12.'*xn2(n+M:-1:n-M);
        y2 = bn21.'*xn1(n+M:-1:n-M)+bn22.'*xn2(n+M:-1:n-M);
        yn1(y_idx) = y1;
        yn2(y_idx) = y2;
        y_idx = y_idx + 1;
        if y_idx > length(yn1)
            break;
        end
        if Update_Tap_or_Not
            switch SE
                case 5
                    en1 = rx_32QAM_Decision(y1) - y1;
                    en2 = rx_32QAM_Decision(y2) - y2;
                case 4
                    en1 = rx_16QAM_Decision(y1) - y1;
                    en2 = rx_16QAM_Decision(y2) - y2;
                case 3
                    en1 = tx_8QAM_mod(rx_8QAM_Decode(y1)) - y1;
                    en2 = tx_8QAM_mod(rx_8QAM_Decode(y2)) - y2;
                case 2
                    en1 = rx_QPSK_Decision(y1) - y1;
                    en2 = rx_QPSK_Decision(y2) - y2;
            end
            if LMS_Plot_Conv_or_Not
                en_list = [en_list,en1];
                bn_list = [bn_list,bn11(M-3)]; 
            end
            bn11 = bn11 + mu*en1*conj(xn1(n+M:-1:n-M));
            bn12 = bn12 + mu*en1*conj(xn2(n+M:-1:n-M));
            bn21 = bn21 + mu*en2*conj(xn1(n+M:-1:n-M));
            bn22 = bn22 + mu*en2*conj(xn2(n+M:-1:n-M));
        end
    end

    if LMS_Plot_Conv_or_Not && Update_Tap_or_Not
        figure; plot(abs(en_list));
        figure; plot(abs(bn_list)); hold on; plot(angle(bn_list));
    end

    yn1 = yn1.';
    yn2 = yn2.';
    
end
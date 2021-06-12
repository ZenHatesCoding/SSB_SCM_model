%% Based on QAMx_LMS_Train
% Update by ZP 11/01/2020
% MIMO
% bnij : input j, output i
%%
% Update by ZP 05/30/2018
% modify to be similar to DP_LMS_Train
%%
% By ZP
% xn, received bit sequence,
% dn, desired bit sequence
% M number of taps
% mu update step
% Sps_in # Sample per symbol of input waveform 
function [bn11,bn12,bn21,bn22] = QAMx_LMS_Train_MIMO(xn1,xn2, dn1,dn2, M, mu, Sps_in,LMS_Plot_Conv_or_Not)
    if size(xn1,2)>1
        xn1 = xn1.';
    end
    
    if size(dn1,2)>1
        dn1 = dn1.';
    end
    
    if size(xn2,2)>1
        xn2 = xn2.';
    end
    
    if size(dn2,2)>1
        dn2 = dn2.';
    end
    
    xn1 = [zeros(M,1);xn1];
    xn2 = [zeros(M,1);xn2];
    dn1 = [zeros(M,1);dn1];
    dn2 = [zeros(M,1);dn2];
    
    
    bn11 = zeros(2*M+1,1);bn11(M+1) = 1; 
    bn12 = zeros(2*M+1,1); 
    bn21 = zeros(2*M+1,1); 
    bn22 = zeros(2*M+1,1);bn22(M+1) = 1; 
    
    dn_idx = 1+M;
    if LMS_Plot_Conv_or_Not
        en_track = [];
        bn_list = [];
    end
    for n = (M+1):Sps_in:length(xn1)-M
        yn1 = bn11.'*xn1(n+M:-1:n-M)+bn12.'*xn2(n+M:-1:n-M);
        yn2 = bn21.'*xn1(n+M:-1:n-M)+bn22.'*xn2(n+M:-1:n-M);
        en1 = (dn1(dn_idx) - yn1);
        en2 = (dn2(dn_idx) - yn2);
        
        if LMS_Plot_Conv_or_Not
            en_track = [en_track en1];
        end
        dn_idx = dn_idx + 1;
        if dn_idx > length(dn1)-M
            break ;
        end
        bn11 = bn11 + mu*en1*conj(xn1(n+M:-1:n-M));
        bn12 = bn12 + mu*en1*conj(xn2(n+M:-1:n-M));
        bn21 = bn21 + mu*en2*conj(xn1(n+M:-1:n-M));
        bn22 = bn22 + mu*en2*conj(xn2(n+M:-1:n-M));

        if LMS_Plot_Conv_or_Not
            bn_list = [bn_list, bn11(M+1)];
        end
    end    
    if LMS_Plot_Conv_or_Not
        figure; plot(abs(en_track));
        figure; plot(abs(bn_list));
    end
end
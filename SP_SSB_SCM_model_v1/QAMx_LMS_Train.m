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
function bn = QAMx_LMS_Train(xn, dn, M, mu, Sps_in,LMS_Plot_Conv_or_Not)
    if size(xn,2)>1
        xn = xn.';
    end
    if size(dn,2)>1
        dn = dn.';
    end
    xn = [zeros(M,1);xn];
    dn = [zeros(M,1);dn];
    
    bn = zeros(2*M+1,1);bn(M+1) = 1; 
    dn_idx = 1+M;
    if LMS_Plot_Conv_or_Not
        en_track = [];
        bn_list = [];
    end
    for n = (M+1):Sps_in:length(xn)-M
%          yn = bn.'*xn(n-M:1:n+M);
        yn = bn.'*xn(n+M:-1:n-M);
        en = (dn(dn_idx) - yn);
        if LMS_Plot_Conv_or_Not
            en_track = [en_track en];
        end
        dn_idx = dn_idx + 1;
        if dn_idx > length(dn)-M
            break ;
        end
        bn = bn + mu*en*conj(xn(n+M:-1:n-M));
%         bn = bn + mu*(en)*conj(xn(n-M:1:n+M));
        if LMS_Plot_Conv_or_Not
            bn_list = [bn_list, bn(M+1)];
        end
    end    
    if LMS_Plot_Conv_or_Not
        figure; plot(abs(en_track));
        figure; plot(abs(bn_list));
    end
end
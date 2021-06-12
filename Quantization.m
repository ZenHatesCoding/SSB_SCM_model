% qin is zero mean
% qout is not normalized
% requires pwr_norm in simulation
function qout = Quantization (qin, nbit)
    nlevel = 2^(nbit-1)-1;
    
    qinmax = max(qin);
    
    qout = round(qin./qinmax*nlevel);
    
        
end
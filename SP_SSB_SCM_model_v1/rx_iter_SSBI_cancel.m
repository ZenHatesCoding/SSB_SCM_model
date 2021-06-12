function y = rx_iter_SSBI_cancel(x,num_iter)
    
    y = x;
    alpha = 3e-3;
    for idx = 1:num_iter
        x_ssb = hilbert(y);
        y = x-alpha*x_ssb.*conj(x_ssb);
        
    end
    
    
    
end
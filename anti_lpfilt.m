% Modified by ZP
% s_rate: sampling rate
% bw: LPF bw
% type: 'gaussian'
% ord: LPF order
function X=anti_lpfilt(X,s_rate,bw,type,ord,cc)

if ~exist('cc','var')
   cc=0; 
end




nrows = size(X,1);
ncols = size(X,2);

Xlen=max(ncols,nrows); % length

if mod(Xlen,2) == 0
    f= linspace(0,s_rate/2,Xlen/2);
else
    f= linspace(0,s_rate/2,(Xlen+1)/2);
end
switch(type)
    
    case 'gaussian'
   
        HLP=exp(-0.5 *log(2)*(f / (bw)).^(2*ord)).';

        if mod(Xlen,2) == 0
            HLP = [HLP;flipud(HLP)];
        else
            HLP = [HLP;flipud(HLP(2:end))];
        end
    case 'bessel'
        [a,b] = besself(ord,2*pi*bw);
        HLP = freqs(a,b,2*pi*f).';
        figure; plot(f/1e9,20*log10(HLP./max(HLP))); grid on;
        if mod(Xlen,2) == 0
            HLP = [HLP;flipud(HLP)];
        else
            HLP = [HLP;flipud(HLP(2:end))];
        end
end


if cc
   
    HLP = conj(HLP);
    
end


if ncols>nrows %horizontal

        X=X.'; % transpose to vertical

        transp = 1;
   
    if nrows == 2 % 2 pol

        npols = 2;
        
    else
        
        npols = 1;
        
    end
        
else  % vertical
    
    transp = 0;
    
    if ncols == 2 % 2 pol
        
        npols = 2;
        
    else
        
        npols = 1;
        
    end
    
end

    
    
    if npols == 2 % 2 pol
    
        X=ifft(fft(X)./[HLP,HLP]);

    else
        
        X=ifft(fft(X)./HLP);

    end
        
    if transp==1  % re-transpose to horizontal

        X=X.';
        
    end
    
    
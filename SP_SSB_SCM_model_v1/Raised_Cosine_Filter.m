%% Principle - Commented by ZP
% *******************************************************************
% RC or RRC is an analog domain filter which shapes the spectrum of the
% input waveform to the output waveform
% Y_a(Omega) = RC_a/RRC_a(Omega)*X_a(Omega)
% x_a(t) and y_a(t) are analog waveforms
% *******************************************************************
% RC filter:
%          cos(pi*beta*t/Tsy)     sin(pi*t/Tsy)
% h_a(t) = ------------------ x -------------  
%          1-(2*beta*t/Tsy)^2       pi*t/Tsy
% 1/Tsy symbol rate; 1/Ts sample rate; Tsy/Ts = sps #samples per symbol
% Notice that h_a(i*Tsy) = 0, i belongs to integer and i!=0
% *******************************************************************
% RRC filter
%             1       sin(pi*t/Tsy*(1-beta))+4*beta*t/Tsy*cos(pi*t/Tsy(1+beta))
% h_a(t) = -------- x ----------------------------------------------------
%          sqrt(Tsy)            pi*t/Tsy*(1-(4*beta*t/Tsy)^2)
% sign doesn't matter using 2 rrc to achieve a rc, -1 will be cancelled
% in this function actually there is a factor -1 comparing to the above
% formula
% |H_rrc| = sqrt(|H_rc|)
% *******************************************************************
% Model: filtering in the Fs sampling system
% Xsy(omega)->D/A->X_a(Omega)->A/D->X(omega)->H(omgea)->Y(omega)->D/A->Y_a(Omega)
%              |                |                                  |
%              Tsy              Ts                                 Ts
% xsy symbol sequence, x,y sample sequence at Fs sampling rate
% y[n] = conv(h[n],x[n]), h[n] = h_a(n*Ts);
% 
% The output is truncated digital filter h, but non-causal,
% needs buffer in real applications 
% 
%% Function
function filter_coef = Raised_Cosine_Filter(Param)

    Time = (-Param.NbSymbols/2:1/Param.OverSamp:Param.NbSymbols/2) + eps;  
    % Time vector
    % number of taps  = Param.OverSamp*(Param.NbSymbols+1) 
    % considering Parm.NbSymbols+1 Symbols, then truncate 
    
    % analog time vector sampled at Fs: t =
    % (-Param.NbSymbos/2*Param.OverSamp:Param.NbSymbols/2*Param.OverSamp)*Ts
    % digtial time vector: n = t/Ts = 
    % -Param.NbSymbos/2*Param.OverSamp:Param.NbSymbols/2*Param.OverSamp
    % Time = n/Param.OverSamp = n/sps = n*Ts/Tsy = t/Tsy
    
    
% %     Freq = (0:Param.NbSamples*Param.FFTPad/2-1)*(Param.OverSamp/2)/(Param.NbSamples*Param.FFTPad/2-1) ;
    
    % Only For plotting purpose
    % maximum is Fs/2 (Param.OverSamp/2, nyquist freq), Dfreq=Freq/(NSamp-1)
    % normalized to Baud rate.
  

    T = 1; % 1 sample per symbol
    a = pi*Time/T;
    
    
    if Param.Type == 1
    %% Raised cosine filter
        RC = sin(a)./a.*(cos(Param.Alpha*.999999*a)./(1-4*(Param.Alpha*.999999)^2 * Time.^2/T^2));
        filter_coef = RC / sum(RC); % normalization
    %% Plot RC filter
%       H_RC = abs(fft(RC,Param.NbSamples*Param.FFTPad));
%       H_RC = H_RC(1:length(H_RC)/2);
%       figure;
%       plot(Freq,20*log10(H_RC));
    elseif Param.Type == 0
    %% Root raised cosine filter

        RRC = 4*Param.Alpha*( cos((1+Param.Alpha*.999999)*a) +...
            sin((1-Param.Alpha*.999999)*a)./(4*Param.Alpha*.999999*Time/T))...
            ./ (pi*T^0.5*((4*Param.Alpha*.999999*Time/T).^2-1));
        filter_coef = RRC / sum(RRC);
    %% PLot RRC filter
%       H_RRC = abs(fft(RRC,Param.NbSamples*Param.FFTPad));
%       H_RRC = H_RRC(1:length(H_RRC)/2);
%     
%       figure
%       plot(RC,'.-'), hold on
%       plot(RRC,'.-r')
% 
%       hold on; plot(Freq,20*log10(H_RRC),'r')  
    end
  
end


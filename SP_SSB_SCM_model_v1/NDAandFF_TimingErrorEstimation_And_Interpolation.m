function [interpolator_Output All_mu_k All_m_k] = NDAandFF_TimingErrorEstimation_And_Interpolation(SignalIn, ProcessingBlockLengthInSamples, Param)
% SignalIn: your signal at 2 Sps;
% ProcessingBlockLengthInSamples: length of blocks in samples
% Param.rollOff: roll-off factor
% Param.originalLengthBeforeImposingSFO: length of SignalIn at 2 Sps

rollOff = Param.rollOff;

i=1;
ProcessingBlock = ((i-1)*ProcessingBlockLengthInSamples+1 : i*ProcessingBlockLengthInSamples);
n = transpose(0:(length(SignalIn(ProcessingBlock))-1));

g = 2*sin(pi*rollOff/2) / (pi*(4-rollOff^2)) / (rollOff/4);
env1 = abs(SignalIn(ProcessingBlock)).^2;
env2 = real(conj(SignalIn(ProcessingBlock)) .* circshift(SignalIn(ProcessingBlock),-1));
estimatedTimingOffset_old = 2 * 1/(2*pi) * angle( g * sum(env1 .* exp(-1i*n*pi)) + sum(env2(1:end-1) .* exp(-1i*pi*(n(1:end-1)-0.5))) );   % relative to the sampling interval

N = length(SignalIn);
f = [(0:N/2-1) (-N/2:-1)].'./(N/2)*0.5;
SignalIn = real(ifft(fft(SignalIn).*exp(1i*2*pi*f*estimatedTimingOffset_old)));
estimatedTimingOffset_old = 0;
% interpolation to correct SFO

% initialization
m_k = ProcessingBlockLengthInSamples + 2;   % basepoint of the sample set
mu_k = estimatedTimingOffset_old;  % fractional delay in units if T_s (sampling interval)
k = ProcessingBlockLengthInSamples + 1;
All_m_k = zeros(1,Param.originalLengthBeforeImposingSFO-10);
All_mu_k = zeros(1,Param.originalLengthBeforeImposingSFO-10);
interpolator_Output = zeros(1,Param.originalLengthBeforeImposingSFO-10);
interpolator_Output(ProcessingBlock) = SignalIn(ProcessingBlock);

All_m_k(ProcessingBlock) = 1:ProcessingBlockLengthInSamples;
All_mu_k(ProcessingBlock) = mu_k;

while m_k <= (length(SignalIn)- 2)  %for k = 1:length(receivedDecimatedSamplesWthSampPhaseErr)
        
% 1) linear interpolation

%     interpolator_Output(k) = SignalIn(m_k) + ...
%         mu_k * (SignalIn(m_k+1) - SignalIn(m_k));
    
% 2) Cubic interpolation  

%     interpolator_Output(k) = SignalIn(m_k+2) * (1/6*mu_k^3-1/6*mu_k) + ...
%         SignalIn(m_k+1) * (-1/2*mu_k^3+1/2*mu_k^2+mu_k) + ...
%         SignalIn(m_k) * (1/2*mu_k^3-mu_k^2-1/2*mu_k+1) + ...
%         SignalIn(m_k-1) * (-1/6*mu_k^3+1/2*mu_k^2-1/3*mu_k);


% 3) piecewise parabolic   (better than qubic)
   
    interpolator_Output(k) = SignalIn(m_k+2) * (1/2*mu_k^2-1/2*mu_k) + ...
        SignalIn(m_k+1) * (-1/2*mu_k^2+3/2*mu_k) + ...
        SignalIn(m_k) * (-1/2*mu_k^2-1/2*mu_k+1) + ...
        SignalIn(m_k-1) * (1/2*mu_k^2-1/2*mu_k);
    
    % check if next block or current block
    
    if rem(m_k,ProcessingBlockLengthInSamples) ~= 0
        m_k = m_k + 1;
    else
        i = i+1;
        if (i*ProcessingBlockLengthInSamples) > length(SignalIn)

            ProcessingBlock = ((i-1)*ProcessingBlockLengthInSamples+1 : length(SignalIn));
            n = transpose(0:(length(SignalIn(ProcessingBlock))-1));
            env1 = abs(SignalIn(ProcessingBlock)).^2;
            env2 = real(conj(SignalIn(ProcessingBlock)) .* circshift(SignalIn(ProcessingBlock),-1));
            estimatedTimingOffset_new = 2 * 1/(2*pi) * angle( g * sum(env1 .* exp(-1i*n*pi)) + sum(env2(1:end-1) .* exp(-1i*pi*(n(1:end-1)-0.5))) );
            
        else
            
            ProcessingBlock = ((i-1)*ProcessingBlockLengthInSamples+1 : i*ProcessingBlockLengthInSamples);
            env1 = abs(SignalIn(ProcessingBlock)).^2;
            env2 = real(conj(SignalIn(ProcessingBlock)) .* circshift(SignalIn(ProcessingBlock),-1));
            estimatedTimingOffset_new = 2 * 1/(2*pi) * angle( g * sum(env1 .* exp(-1i*n*pi)) + sum(env2(1:end-1) .* exp(-1i*pi*(n(1:end-1)-0.5))) );
            
        end
        Adddrop= floor(mu_k + 1 + SAW(estimatedTimingOffset_new - estimatedTimingOffset_old));
        m_k = m_k + Adddrop;
        mu_k = mod(mu_k + 1 + SAW(estimatedTimingOffset_new - estimatedTimingOffset_old) , 1);   % new mu_k
%  
%         AddDrop(jj) = ceil( (SFO_Vector(jj-1)-dt+1)/2);
%         dt = dt+ceil( (SFO_Vector(jj-1)-dt+1)/2)-1;
% m_k = m_k + Adddrop;
%          SFO_Vector(jj)= dt;
%         
         estimatedTimingOffset_old = estimatedTimingOffset_new;

    end
     
    All_m_k(k) = m_k;
    All_mu_k(k) = mu_k;
    k = k + 1;
end

interpolator_Output = interpolator_Output.';


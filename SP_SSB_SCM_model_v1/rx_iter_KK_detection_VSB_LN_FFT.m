% ***************************************************
% M: number of taps used in overlap and save algorithm
% 1 iteration only
% ***************************************************
function Rx_Time_Data = rx_iter_KK_detection_VSB_LN_FFT(Rx_I,ParamRxDSP,...
    DTime,ParamVSB, ParamFib, ParamSig, ParamControl)
    
    %% prepare parameters
    totalBaudRate = sum(ParamSig.SC_Baud_Rate);
    roll_off = ParamSig.roll_off;
    GuardBand = ParamSig.GuardBand;
    
    Fiber_Length = ParamFib.FiberLength; 
    Beta2 = ParamFib.Beta2_ref;
    
    %% get VFreq
    
    Nfft = round(ParamRxDSP.Nfft*ParamRxDSP.KKoverSamp*28/56);
%     Nfft = ParamRxDSP.Nfft;
    M = ParamRxDSP.FAiterKK_tap;
    L = Nfft-M;
    num_block = 2*ceil(length(Rx_I)/L/2);
%     num_block = ceil(length(Rx_I)/L);
    Rx_I = [Rx_I, zeros(1,num_block*L-length(Rx_I))];
    
    
    
    MaxFreq = 0.5/DTime;
    K = 16;
    VFreq = linspace(-MaxFreq,MaxFreq,Nfft*K);
    VFreq = VFreq(1:K:end);
    VOmeg = 2*pi*VFreq;	
    
    Disper = -(1j/2)*Beta2*VOmeg.^2*Fiber_Length; 

    Freq_min = GuardBand;
    Freq_max = GuardBand+(1+roll_off)*totalBaudRate;


%%                  
    switch ParamControl.OBPF_option
        case 1
            Attenuation_dB = zeros(size(VFreq)); 
            Slope_Stop = ParamVSB.Opt_Flt_offset; 
            Slope_Start =Slope_Stop-ParamVSB.Opt_Flt_Suppression_dB/ParamVSB.Opt_Flt_Slope_dBper10GHz*10e9;
            if min(VFreq) < Slope_Start
                Attenuation_dB(VFreq<=Slope_Start)= -ParamVSB.Opt_Flt_Suppression_dB;
                length_slope = length(Attenuation_dB(VFreq>Slope_Start & VFreq <=Slope_Stop));
                Attenuation_dB(VFreq>Slope_Start & VFreq <=Slope_Stop) = linspace(-ParamVSB.Opt_Flt_Suppression_dB,0,length_slope);
            else
                ParamVSB.Opt_Flt_Suppression_dB = (Slope_Stop-min(VFreq))/10e9*ParamVSB.Opt_Flt_Slope_dBper10GHz;
                length_slope = length(Attenuation_dB(VFreq <=Slope_Stop));
                Attenuation_dB(VFreq <=Slope_Stop) = linspace(-ParamVSB.Opt_Flt_Suppression_dB,0,length_slope);
            end
            
            if ParamControl.LSB_or_RSB
                Attenuation_dB = fliplr(Attenuation_dB);
            end
        case 5
            Mat = csvread('LumentumFilter.csv');
            WL = Mat(:,1)*1e-9;
            Pwr = -Mat(:,2);
            c = 3e8;
            f = c./WL;
            f_center = c/1541.9e-9;

            Freq_Vector_new = VFreq+f_center+ParamVSB.Opt_Flt_offset;
            Pwr_new = interp1(f,Pwr,Freq_Vector_new); 
            Attenuation_dB = Pwr_new-max(Pwr_new);
            if ParamControl.LSB_or_RSB == 0
                Attenuation_dB = fliplr(Attenuation_dB);
            end

            Attenuation_dB(isnan(Attenuation_dB)) =min(Attenuation_dB);
    end

    attenuation = 10.^(0.05*Attenuation_dB);

    H2 = fliplr(attenuation);

    

    H2 = H2.*exp(-1*Disper);
    H2(VFreq<=0) = 1;

    H1 = attenuation;
    H1 = H1.*exp(1*Disper);
    H1(VFreq<=0) = 1;
    
    
    bn = real(sqrt(Rx_I));
    bn = [zeros(1,M),bn];

%     figure;
%     plot(10*log10(abs(H1.*(VFreq>Freq_min &VFreq <=Freq_max)))); hold on;
%     plot(10*log10(abs(H2.*(VFreq>Freq_min &VFreq <=Freq_max)))); hold on;
%     plot(10*log10(abs((H1+H2).*(VFreq>Freq_min &VFreq <=Freq_max))));

   
    if ParamControl.FFT_sharing_or_Not == 0
        for idx = 1:ParamRxDSP.numKKiter
            imS1 = [];
            imS2 = [];
            for i = 1:num_block
                Input_Samples_block = bn((i-1)*L+1:(i-1)*L+Nfft);

                S1_F = fftshift(fft(Input_Samples_block));

%                 alpha = 2; % this is magic ! I don't really understand this.
%                 S1_F = alpha.*S1_F./(H2+H1);
%                 S1_F = S1_F.*(VFreq>Freq_min &VFreq <=Freq_max);
    
                F1 = 2*(VFreq>Freq_min &VFreq <=Freq_max).*H1./(H2+H1);
                F2 = 2*(VFreq>Freq_min &VFreq <=Freq_max).*H2./(H2+H1);

                imS1_block = imag(ifft(ifftshift(S1_F.*F1)));
                imS2_block = imag(ifft(ifftshift(S1_F.*F2)));
                
                

                
                imS1  = [imS1,imS1_block(M+1:end)];
                imS2  = [imS2,imS2_block(M+1:end)];
            end
            bn = real(sqrt(Rx_I-(imS1-imS2).^2));

            bn = [zeros(1,M),bn];

        end
       
        
    else        
        for idx = 1:ParamRxDSP.numKKiter
            imS1 = [];
            imS2 = [];
            diff_im = [];
            
            for i = 1:2:num_block
                Input_Samples_block = bn((i-1)*L+1:(i-1)*L+Nfft);
                
                if i == num_block
                    Input_Samples_block2 = zeros(1,Nfft);
                else
                    Input_Samples_block2 = bn((i)*L+1:(i)*L+Nfft);
                
                end
                
                Freq_Samples_block = fftshift(fft(Input_Samples_block+1i*Input_Samples_block2));
                
                Freq_Samples_block_cr = real(Freq_Samples_block);
                Freq_Samples_block_ci = imag(Freq_Samples_block);
                
                Freq_Samples_block_r = 1/2*(Freq_Samples_block_cr + fliplr(Freq_Samples_block_cr)); 
                Freq_Samples_block_i = 1/2*(Freq_Samples_block_ci - fliplr(Freq_Samples_block_ci)); 

                X_block1 =  Freq_Samples_block_r + 1i*Freq_Samples_block_i;

                Freq_Samples_block_r = 1/2*(Freq_Samples_block_ci + fliplr(Freq_Samples_block_ci)); 
                Freq_Samples_block_i = -1/2*(Freq_Samples_block_cr - fliplr(Freq_Samples_block_cr)); 

                X_block2 = Freq_Samples_block_r + 1i*Freq_Samples_block_i;
                
               
                
                
                F1 = (VFreq>Freq_min &VFreq <=Freq_max).*H1./(H2+H1);
                F2 = (VFreq>Freq_min &VFreq <=Freq_max).*H2./(H2+H1);
                
                F11 = -1i*F1+1i*conj(fliplr(F1));
                F12 = 1*F1-1*conj(fliplr(F1));
                F21 = -1i*F2+1i*conj(fliplr(F2));
                F22 = 1*F2-1*conj(fliplr(F2));

%                 Output_Freq_Samples_block1 = (X_block1.*F11+X_block2.*F12);
%                 Output_Freq_Samples_block2 = (X_block1.*F21+X_block2.*F22);
                
                
%                 Output_Samples_block1 = ifft(ifftshift(Output_Freq_Samples_block1));
%                 Output_Samples_block2 = ifft(ifftshift(Output_Freq_Samples_block2));
%                 
%                 imS1  = [imS1,real(Output_Samples_block1(M+1:end)),imag(Output_Samples_block1(M+1:end))];
%                 imS2  = [imS2,real(Output_Samples_block2(M+1:end)),imag(Output_Samples_block2(M+1:end))];
%               
                Output_Freq_Samples_block = X_block1.*(F11-F21)+X_block2.*(F12-F22);
%                 Output_Freq_Samples_block = Output_Freq_Samples_block1-Output_Freq_Samples_block2;
                Output_Samples_block = ifft(ifftshift(Output_Freq_Samples_block));
                diff_im  = [diff_im,real(Output_Samples_block(M+1:end)),imag(Output_Samples_block(M+1:end))];


            end
            
            bn = real(sqrt(Rx_I-(diff_im).^2));

            bn = [zeros(1,M),bn];

        end 
    
        
    end
    
    Rx_Time_Data = bn(M+1:end);
end 
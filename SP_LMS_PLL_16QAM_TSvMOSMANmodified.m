% -------------------------------------------------------------------------
% PMD post-compensation
% By : Benoit Chatelain
% July 10 2009, McGill University

% v1
% -------------------------------------------------------------------------

function SignalXOut = SP_LMS_PLL_16QAM_TSvMOSMANmodified(Param,SignalX,Training_Symbols_X)
    [numRow, numColumn] = size(SignalX);
        if numColumn == 1
            SignalX = SignalX.';
        end 
    TS_Block_Length = length(Training_Symbols_X);

    g = Param.PLLg;
    PLL_Error_X = zeros(1, round(length(SignalX)/2));
    PLL_Phase_X = PLL_Error_X;

    mu1 = Param.Trainmu;
    mu2 = Param.DDmu;

    hxx = zeros(1,Param.NbTapsLMS) + 1i*zeros(1,Param.NbTapsLMS);

        hxx((Param.NbTapsLMS+1)/2) = 1;


    AllError  = zeros(round(length(SignalX)/2),1);
    AllCoefs = zeros(round(length(SignalX)/2),Param.NbTapsLMS);

    SignalXOut = zeros(1,round(length(SignalX)/2));

    SignalXPadded = [zeros(1,(Param.NbTapsLMS-1)/2) SignalX zeros(1,(Param.NbTapsLMS-1)/2)];

    currentIndex = 1;

    for aa = Param.NbTapsLMS:2:length(SignalXPadded) % Equalizer taps setting


            % Align input reference signal
            xx = SignalXPadded(aa-Param.NbTapsLMS+1:aa);

            xx = fliplr(xx);

            % Compute output signal for each polarization
            xout = hxx * xx.';

            SignalXOut(currentIndex) = xout * exp(-1j*PLL_Phase_X(currentIndex));

           % Decision

            if currentIndex <= TS_Block_Length


                xout_CPR_Decision = Training_Symbols_X(currentIndex);

                % Compute error signal for each polarization

                Ex = xout_CPR_Decision * exp(1j*PLL_Phase_X(currentIndex)) - xout;

    %             Ex = xout_CPR_Decision - SignalXOut(currentIndex);
                % Update PMD filter coefficients (LMS algorithm)
                hxx = hxx + mu1*Ex*conj(xx);


            else



                xout_CPR_Decision = rx_16QAM_Decision(SignalXOut(currentIndex));
    %             xout_CPR_Decision = rx_QPSK_Decision(SignalXOut(currentIndex));
                Ex = xout_CPR_Decision * exp(1j*PLL_Phase_X(currentIndex)) - xout;
    %             Ex = xout_CPR_Decision - SignalXOut(currentIndex);

                % Update PMD filter coefficients (LMS algorithm)
                hxx = hxx + mu2*Ex*conj(xx);

            end


            %PLL
    %         PLL_Error_X(currentIndex) = angle(SignalXOut(currentIndex).*conj(rx_16QAM_Decision(SignalXOut(currentIndex))));
    %         PLL_Error_Y(currentIndex) = angle(SignalYOut(currentIndex).*conj(rx_16QAM_Decision(SignalYOut(currentIndex))));

            PLL_Error_X(currentIndex) = angle(SignalXOut(currentIndex).*conj(xout_CPR_Decision));


            PLL_Phase_X(currentIndex+1) = g*PLL_Error_X(currentIndex)+PLL_Phase_X(currentIndex);




            AllError(currentIndex,1) = abs(Ex)^2;

            AllCoefs(currentIndex,:) =  hxx;


            currentIndex = currentIndex +1;
    end

    if Param.DisplayFigures == 1

        figure
        plot(real(AllCoefs(:,(Param.NbTapsLMS+1)/2,1)),'r-');
        legend('hxx');
    end

    if Param.DisplayFigures == 1
        figure
        plot(real(hxx)), hold on
        plot(imag(hxx),'r')
        title('hxx')

    end

    if numColumn == 1
        SignalXOut = SignalXOut.';
    end
end
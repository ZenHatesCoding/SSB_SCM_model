
 
clc; 
clear;
close all;

%% Initialization
SP_VSB_Simulator_Control_Param_Init;

ParamControl.CSPR_tuning_case = 1; % 1 Linear model 2 Parallel Path 3 bias tuning
SP_VSB_Simulator_System_Param_Init;


Path = 'D:\xzp\Project\VSB_filter_requirement\';


numtrial = 10;
formatspec = [];
for i = 1:numtrial
    formatspec = [formatspec,'%12.4e '];
end

formatspec = [formatspec, '\n'];

for KK_option = [3]
    ParamControl.KK_option = KK_option; %0 normal KK, 1 KK w/o upsamp
    ParamControl.CD_Compensation_or_Not = 1;
    for KK_iter = [1]
        ParamRxDSP.numKKiter = KK_iter;

        if ParamControl.KK_option == 3
            ParamControl.CD_Compensation_or_Not = 1;
            TxT_name0 = [Path,'Drift_study_CSPRswp_FAitKK_PAM_VSB_iter_',num2str(KK_iter),'_'];
        else
            ParamControl.CD_Compensation_or_Not = 1;
            TxT_name0 = [Path,'cKK_PAM_VSB_'];
        end

        for Length = [40e3]
            ParamFib.FiberLength = Length;
            zp_Fib_dispersion_Param_Init;
            TxT_name = [TxT_name0,num2str(Length/1e3),'km.txt'];
            screen_message = ['Length_',num2str(Length/1e3),'km_'];
            for slope = [7,11,15]
                ParamVSB.Opt_Flt_Slope_dBper10GHz = slope;
                screen_message1 = [screen_message, 'slope_',num2str(slope),'dBper10GHz_'];
                for OSNR = [32]
                    ParamChan.SNR = OSNR  + ...
                        10*log10(ParamGen.NoiseReferenceBW/ParamSig.Baud_Rate);
                    screen_message2 = [screen_message1,'OSNR_',num2str(OSNR),'dB_'];
                    for CSPR = [6,7,8,9]
                        ParamSys.CSPR_dB = CSPR;
                        screen_message3 = [screen_message2,'CSPR_',num2str(CSPR),'dB'];
                        for wavelength_drift = -2:0.2:2
                            ParamVSB.Opt_Flt_drift = wavelength_drift*1e9;
                            screen_message4 = [screen_message3,'WL_drift_',num2str(wavelength_drift),'GHz'];
                            BER_list = zeros(1,numtrial);
                            for i = 1:numtrial
                                SP_SSB_PhysicalLayerSimulator;
                                disp('CSPR:')
                                disp(ParamPhysicalModel.Measured_CSPR);
                                disp('BER:')
                                disp(ParamPhysicalModel.BER_avg);
                                BER_list(i) = ParamPhysicalModel.BER_avg;
                            end

                            fileID = fopen(TxT_name,'at');

                            fprintf(fileID,formatspec,BER_list);
                            fprintf(fileID,[screen_message4,'\n']);
                            fclose(fileID);
                            disp(screen_message4);
                        end
                    end
                end
            end
        end
    end
end




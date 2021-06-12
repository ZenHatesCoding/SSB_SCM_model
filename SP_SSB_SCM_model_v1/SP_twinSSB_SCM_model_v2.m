
 
clc; 
clear;
close all;

%% Initialization
SP_twinSSB_Simulator_Control_Param_Init;

%---------
% disable all plots during sweep
ParamControl.Plot_opt_filter_or_Not = 0;
ParamControl.Plot_Spectrum_or_Not = 0;
ParamControl.Plot_CrossCorr_or_Not = 0;
ParamControl.Equalization_Plot_Conv_or_Not = 0;
ParamControl.Plot_Constellation_or_Not = 0;
ParamControl.Plot_Error_or_Not = 0;
ParamControl.Plot_Power_Pie_or_Not = 0;
%---------
SP_twinSSB_Simulator_System_Param_Init;


Path = 'D:\xzp\Project\twinSSB_filter_requirement\';


numtrial = 10;
formatspec = [];
for i = 1:numtrial
    formatspec = [formatspec,'%12.4e '];
end

formatspec = [formatspec, '\n'];

for KK_option = [3 0]
    ParamControl.KK_option = KK_option; %0 normal KK, 1 KK w/o upsamp
    ParamControl.CD_Compensation_or_Not = 1;
    if ParamControl.KK_option == 3
        ParamControl.CD_Compensation_or_Not = 0;
        TxT_name0 = [Path,'BothFilter_FAitKK_twinSB_'];
    else
        ParamControl.CD_Compensation_or_Not = 1;
        TxT_name0 = [Path,'BothFilter_cKK_twinVSB_'];
    end

    for Length = [40e3]
        ParamFib.FiberLength = Length;
        zp_Fib_dispersion_Param_Init;
        TxT_name = [TxT_name0,num2str(Length/1e3),'km.txt'];
        screen_message = ['Length_',num2str(Length/1e3),'km_'];
        for slope = 19:-4:3
            ParamVSB.Opt_Flt_Slope_dBper10GHz1 = slope;
            ParamVSB.Opt_Flt_Slope_dBper10GHz2 = slope;
            screen_message1 = [screen_message, 'slope_',num2str(slope),'dBper10GHz_'];
            for OSNR = 40:-1:26
                ParamChan.SNR = OSNR  + ...
                    10*log10(ParamGen.NoiseReferenceBW/ParamSig.Baud_Rate);
                screen_message2 = [screen_message1,'OSNR_',num2str(OSNR),'dB_'];
                for CSPR = 6:14
                    ParamSys.CSPR_dB = CSPR;
                    screen_message3 = [screen_message2,'CSPR_',num2str(CSPR),'dB'];
                    BER_list = zeros(1,numtrial);
                    for i = 1:numtrial
                        SP_twinSSB_PhysicalLayerSimulator;
                        disp('CSPR:')
                        disp(ParamPhysicalModel.Measured_CSPR);
                        disp('BER:')
                        disp(ParamPhysicalModel.BER_avg);
                        BER_list(i) = ParamPhysicalModel.BER_avg;
                    end

                    fileID = fopen(TxT_name,'at');

                    fprintf(fileID,formatspec,BER_list);
                    fprintf(fileID,[screen_message3,'\n']);
                    fclose(fileID);
                    disp(screen_message3);
                end
            end
        end
    end
end




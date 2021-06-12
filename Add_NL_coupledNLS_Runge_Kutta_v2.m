%% doesn't work, maybe the step size is not small enough ?

% dA1/dZ = f(A1,z) 
%        = -alpha/2*A1 - 1i*Beta0*A1- 1i*gamma*[(|A1|^2+2*sum(|Ai|^2))A1 + ...
%          A2*A3*conj(A4)]
% deltaK: phase mismatch


% Runge Kutta method
% dy/dz = f(z,y), y0 z0 known
% k1 = f(z0,y0)
% k2 = f(z0+dZ/2, y0 + dZ*k1/2)
% k3 = f(z0+dZ/2, y0 + dZ*k2/2)
% k4 = f(z0+dZ, y0 + dZ*k3)

% z1 = z0+dZ
% y1 = y0 + 1/6*dZ*(k1+2*k2+2*k3+k4)

function [Tx_Time_Data_CWDM_out,Z1] = Add_NL_coupledNLS_Runge_Kutta_v2(Tx_Time_Data_CWDM,ParamFib,dZ,Z0) 
    Y0 = Tx_Time_Data_CWDM;

    k11 = Func1(Y0{1},Y0{2},Y0{3},Y0{4},ParamFib,Z0);
    k12 = Func2(Y0{1},Y0{2},Y0{3},Y0{4},ParamFib,Z0);
    k13 = Func3(Y0{1},Y0{2},Y0{3},Y0{4},ParamFib,Z0);
    k14 = Func4(Y0{1},Y0{2},Y0{3},Y0{4},ParamFib,Z0);
    
    
    k21 = Func1(Y0{1} + dZ/2*k11, Y0{2} + dZ/2*k12, Y0{3} + dZ/2*k13, Y0{4} + dZ/2*k14,...
                ParamFib,Z0 + dZ/2);
    k22 = Func2(Y0{1} + dZ/2*k11, Y0{2} + dZ/2*k12, Y0{3} + dZ/2*k13, Y0{4} + dZ/2*k14,...
                ParamFib,Z0 + dZ/2);
    k23 = Func3(Y0{1} + dZ/2*k11, Y0{2} + dZ/2*k12, Y0{3} + dZ/2*k13, Y0{4} + dZ/2*k14,...
                ParamFib,Z0 + dZ/2);            
    k24 = Func4(Y0{1} + dZ/2*k11, Y0{2} + dZ/2*k12, Y0{3} + dZ/2*k13, Y0{4} + dZ/2*k14,...
                ParamFib,Z0 + dZ/2);
            
    k31 = Func1(Y0{1} + dZ/2*k21, Y0{2} + dZ/2*k22, Y0{3} + dZ/2*k23, Y0{4} + dZ/2*k24,...
                ParamFib,Z0 + dZ/2);
    k32 = Func2(Y0{1} + dZ/2*k21, Y0{2} + dZ/2*k22, Y0{3} + dZ/2*k23, Y0{4} + dZ/2*k24,...
                ParamFib,Z0 + dZ/2);
    k33 = Func3(Y0{1} + dZ/2*k21, Y0{2} + dZ/2*k22, Y0{3} + dZ/2*k23, Y0{4} + dZ/2*k24,...
                ParamFib,Z0 + dZ/2);
    k34 = Func4(Y0{1} + dZ/2*k21, Y0{2} + dZ/2*k22, Y0{3} + dZ/2*k23, Y0{4} + dZ/2*k24,...
                ParamFib,Z0 + dZ/2);
            
    k41 = Func1(Y0{1} + dZ*k31, Y0{2} + dZ*k32, Y0{3} + dZ*k33, Y0{4} + dZ*k34,...
                ParamFib,Z0 + dZ);
    k42 = Func2(Y0{1} + dZ*k31, Y0{2} + dZ*k32, Y0{3} + dZ*k33, Y0{4} + dZ*k34,...
                ParamFib,Z0 + dZ);
    k43 = Func3(Y0{1} + dZ*k31, Y0{2} + dZ*k32, Y0{3} + dZ*k33, Y0{4} + dZ*k34,...
                ParamFib,Z0 + dZ);
    k44 = Func4(Y0{1} + dZ*k31, Y0{2} + dZ*k32, Y0{3} + dZ*k33, Y0{4} + dZ*k34,...
                ParamFib,Z0 + dZ);
            
    Z1 = Z0 + dZ;
    Y1{1} = Y0{1} + 1/6*dZ*(k11+2*k21+2*k31+k41);
    Y1{2} = Y0{2} + 1/6*dZ*(k12+2*k22+2*k32+k42);
    Y1{3} = Y0{3} + 1/6*dZ*(k13+2*k23+2*k33+k43);
    Y1{4} = Y0{4} + 1/6*dZ*(k14+2*k24+2*k34+k44);
    
    Tx_Time_Data_CWDM_out = Y1;
            
end

function F1 = Func1(A1,A2,A3,A4,ParamFib,Z)
    alpha = ParamFib.Loss_alpha;
    gamma = ParamFib.gamma;
    FWM = ParamFib.FWM_Enable;
    dFWM = ParamFib.dFWM_Enable;

    F1 = -alpha/2*A1 ...
         -1i*gamma*(...
         (abs(A1).^2+...
         2*(abs(A2).^2+abs(A3).^2 + abs(A4).^2)).*A1 + ...
         FWM*2*A2.*A3.*conj(A4)+ ...
         dFWM*1*A2.*A2.*conj(A3));
end
function F2 = Func2(A1,A2,A3,A4,ParamFib,Z)
    alpha = ParamFib.Loss_alpha;
    gamma = ParamFib.gamma;

    FWM = ParamFib.FWM_Enable;
    dFWM = ParamFib.dFWM_Enable;

    F2 = -alpha/2*A2  ...
         -1i*gamma*(...
         (abs(A2).^2+...
         2*(abs(A1).^2+abs(A3).^2 + abs(A4).^2)).*A2 + ...
         FWM*2*A1.*A4.*conj(A3)+ ...
         dFWM*1*A3.*A3.*conj(A4)+...
         dFWM*2*A1.*A3.*conj(A2));
end
function F3 = Func3(A1,A2,A3,A4,ParamFib,Z)
    alpha = ParamFib.Loss_alpha;
    gamma = ParamFib.gamma;
   
    FWM = ParamFib.FWM_Enable;
    dFWM = ParamFib.dFWM_Enable;

    F3 = -alpha/2*A3 ...
         -1i*gamma*(...
         (abs(A3).^2+...
         2*(abs(A1).^2+abs(A2).^2 + abs(A4).^2)).*A3 + ...
         FWM*2*A1.*A4.*conj(A2)+ ...
         dFWM*1*A2.*A2.*conj(A1)+...
         dFWM*2*A2.*A4.*conj(A3));
end
function F4 = Func4(A1,A2,A3,A4,ParamFib,Z)
    alpha = ParamFib.Loss_alpha;
    gamma = ParamFib.gamma;
    
    FWM = ParamFib.FWM_Enable;
    dFWM = ParamFib.dFWM_Enable;

    F4 = -alpha/2*A4 ...
         -1i*gamma*(...
         (abs(A3).^2+...
         2*(abs(A1).^2+abs(A2).^2 + abs(A3).^2)).*A4 + ...
         FWM*2*A2.*A3.*conj(A1)+ ...
         dFWM*1*A3.*A3.*conj(A2));
end
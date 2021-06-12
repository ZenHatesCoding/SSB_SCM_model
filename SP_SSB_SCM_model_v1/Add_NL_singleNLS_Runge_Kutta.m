% dA/dZ = f(A1,z) 
%        = -alpha/2*A - 1i*gamma*|A|^2*A1 + 

% Runge Kutta method
% dy/dz = f(z,y), y0 z0 known
% k1 = f(z0,y0)
% k2 = f(z0+dZ/2, y0 + dZ*k1/2)
% k3 = f(z0+dZ/2, y0 + dZ*k2/2)
% k4 = f(z0+dZ, y0 + dZ*k3)

% z1 = z0+dZ
% y1 = y0 + 1/6*dZ*(k1+2*k2+2*k3+k4)

function [Y1,Z1] = Add_NL_singleNLS_Runge_Kutta(Y0,ParamFib,dZ,Z0) 

    k1 = Func(Y0,ParamFib);
    
    k2 = Func(Y0 + dZ*k1/2,ParamFib);
    
    k3 = Func(Y0 + dZ*k2/2,ParamFib);
    
    k4 = Func(Y0 + dZ*k3,ParamFib); 
    
    Z1 = Z0 + dZ;
    Y1 = Y0 + 1/6*dZ*(k1+2*k2+2*k3+k4);
            
end

function F = Func(A,ParamFib)
    alpha = ParamFib.Loss_alpha;
    gamma = ParamFib.gamma;
    F = -alpha/2*A -1i*gamma*(abs(A).^2).*A;
end

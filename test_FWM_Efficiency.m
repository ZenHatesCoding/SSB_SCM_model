deltaK_nd = ParamFib.Beta0(1)+ParamFib.Beta0(4)-ParamFib.Beta0(3)-ParamFib.Beta0(2);
deltaK_d = ParamFib.Beta0(3)+ParamFib.Beta0(3)-ParamFib.Beta0(4)-ParamFib.Beta0(2);
deltaK_nd2 = ParamFib.Beta0(1)+ParamFib.Beta0(3)-ParamFib.Beta0(2)-ParamFib.Beta0(2);


L = ParamFib.FiberLength;
alpha = ParamFib.Loss_alpha;

phase = (alpha + 1i*deltaK_nd)*L;
eta1 = abs((1-exp(-phase))/phase)^2

phase = (alpha + 1i*deltaK_d)*L;
eta2 = abs((1-exp(-phase))/phase)^2

phase = (alpha + 1i*deltaK_nd2)*L;
eta3 = abs((1-exp(-phase))/phase)^2


phase = alpha*L;
eta0 = abs((1-exp(-phase))/phase)^2

Pwr_diff = 10*log10(eta0./[eta1,eta2,eta3])

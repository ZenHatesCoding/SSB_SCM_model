function phase = Add_PhaseNoise(sizeXin,Sample_Rate,delta_mu)

    phase_sigma_square = 2*pi*delta_mu*1/Sample_Rate;
    phase_noise = sqrt(phase_sigma_square)*randn(sizeXin);

    phase = zeros(sizeXin);
    if sizeXin(1)>sizeXin(2) 
        LengthXin = sizeXin(1);
    else
        LengthXin = sizeXin(2);
    end
    for idx = 2:LengthXin
        phase(idx) = phase(idx-1)+phase_noise(idx);
    end
end


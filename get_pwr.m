function Pwr = get_pwr(x_in)
    Pwr = sum(norm(x_in)^2)/length(x_in);
end
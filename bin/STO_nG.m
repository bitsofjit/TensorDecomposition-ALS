function [STO_1G, STO_3G] = STO_nG(x)

r = norm(x);

fit_coeff_3G = [0.0835; 0.2678; 0.2769];

exp_funcs_3G = [exp(-0.1689*r^2); exp(-0.6239*r^2); exp(-3.4253*r^2)];

STO_1G = 0.3696*exp(-0.4166*r^2);

STO_3G = dot(fit_coeff_3G, exp_funcs_3G);

end

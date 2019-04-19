function y = intial_cond_soliton(x, c, t)
% computes the soliton initial condition of KdV
% phi_sol(.,t)
    %sech_tmp = sech(0.5 * sqrt(c) * (x-c*t));
    y = 0.5*c*sech(0.5*sqrt(c)*(x-c*t)).^2;
end
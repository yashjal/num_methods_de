function y = intial_cond_soliton2(x, c1, c2, t)
% computes the soliton initial condition of KdV
% phi_sol(.,t)
    y = 0.5*c1*sech(0.5*sqrt(c1)*((x-10)-c1*t)).^2+0.5*c2*sech(0.5*sqrt(c2)*((x+10)-c2*t)).^2;
end
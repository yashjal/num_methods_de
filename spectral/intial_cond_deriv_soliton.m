function y = intial_cond_deriv_soliton(x, c)
% computes the derivative of soliton initial condition of KdV
% d_t phi_sol(.,t=0)
    y = 0.5 * (c^(5/2) * sinh(0.5*sqrt(c)*x)) ./ ((cosh(0.5*sqrt(c)*x)).^3);
end
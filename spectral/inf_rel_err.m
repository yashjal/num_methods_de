function y = inf_rel_err(f, phi)
% computes the inf relative error between the interpolant and function over fine grid
    y = norm(f-phi,inf)/norm(f,inf);
end
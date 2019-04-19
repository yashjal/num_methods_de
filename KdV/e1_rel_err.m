function y = e1_rel_err(f, phi, L)
% computes the e1 relative error between the interpolant and function over fine grid
    N = size(f,2);
    h = L/N;
    y = h*norm(f-phi,1)/(h*norm(f,1));
end
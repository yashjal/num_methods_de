function y = e2_rel_err(f, phi, L)
% computes the e2 relative error between the interpolant and function over fine grid
    N = size(f,2);
    h = L/N;
    %y = sqrt(h*(norm(f-phi)^2))/sqrt(h*(norm(f)^2));
    y = sqrt(h)*norm(f-phi)/(sqrt(h)*norm(f));
end

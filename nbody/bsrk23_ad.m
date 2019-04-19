function [t,y,del,tall,xall] = bsrk23_ad(f,tspan,y0,tol)
% adaptive Bogacki-Shampine RK method
    p = 2; % convergence order of lower order estimate
    tn = tspan(1);
    T = tspan(2);
    h = 0.3; % arbitrary initialization of time step
    yn = y0;
    k4 = f(tn,yn);
    t = tn; y = yn';
    tall = tn; xall = yn(1);
    del = h;
    while tn ~= T
        tn1 = tn + h;
        if tn1 > T
            tn1 = T;
        end
        k1 = k4;
        k2 = f(tn+0.5*h, yn+0.5*h*k1);
        k3 = f(tn+0.75*h, yn+0.75*h*k2);
        yn1 = yn + (2/9)*h*k1 + (1/3)*h*k2 + (4/9)*h*k3;
        k4 = f(tn1,yn1);
        z = yn + (7/24)*h*k1 + 0.25*h*k2 + (1/3)*h*k3 + (1/8)*h*k4;
        e0 = h*tol/T; ek = norm(z-yn1);
        if ek < e0
            h = h * min(5, 0.8*(e0/ek)^(1/(p+1)));
        elseif ek > e0
            h = 0.8*h*(e0/ek)^(1/p);
            tall = [tall;tn1]; xall = [xall;yn1(1)];
            continue;
        end
        tn = tn1; yn = yn1;
        y = [y;yn']; t = [t;tn]; del = [del;h]; 
        tall = [tall;tn]; xall = [xall;yn(1)];
    end
end
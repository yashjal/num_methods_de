function [t,y] = bsrk23(f,tspan,y0)
% Bogacki-Shampine RK method
    tn = tspan(1);
    yn = y0;
    k4 = f(tn,yn);
    t = tspan; % assuming tspan consists all the points
    N2 = length(yn);
    y = zeros(length(tspan),N2);
    y(1,:) = yn';
    count = 2;
    for tn1=(tspan(2:end))'
        h = tn1-tn;
        k1 = k4;
        k2 = f(tn+(0.5*h), yn+(0.5*h*k1));
        k3 = f(tn+(0.75*h), yn+(0.75*h*k2));
        yn1 = yn + (2*k1 + 3*k2 + 4*k3) * h/9;
        k4 = f(tn1,yn1);
        %z = yn + (7/24)*h*k1 + 0.25*h*k2 + (1/3)*h*k3 + (1/8)*h*k4;
        tn = tn1;
        yn = yn1;
        y(count,:) = yn';
        count = count + 1;
    end
    
end


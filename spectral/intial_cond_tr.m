function y = intial_cond_tr(x)
% computes the triangle wave initial condition
    n = size(x,2);
    y = zeros(1,n);
    for i=1:n
        if x(i) < 0 && x(i) >= -10
            y(i) = 10 + x(i);
        elseif x(i) >= 0 && x(i) <= 10
            y(i) = 10 - x(i);
        end
    end
end
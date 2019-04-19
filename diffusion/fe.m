function u = fe(A,un,f,dt)
% forward euler step in time
% f = f(xn,tn)
u = un + dt*(A*un+f);
end


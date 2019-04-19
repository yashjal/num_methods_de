function dydt = nbodyodef(t,y,n,G,m,r0)
% rhs f for N-body autonomous ode
    N = length(y)/2; % this N is actually N*n
    dydt = [y((N+1):2*N);zeros(N,1)];
    for i=1:n:N        
        dydt(i+N:n+i+N-1) = F(i,y(1:N),n,G,m,r0);
    end
end


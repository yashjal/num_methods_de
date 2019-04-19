% 1.1.1 hw5 elliptic pde
% solving for u(x,0) and estimating order of accuracy

d = @(x) 2 + cos(pi*x); % d(x)
usol = @(x,t) exp(-pi^2*t/4)*cos(0.5*pi*x); % manufactured soln
fx = @(x) 0.5*pi^2*cos(pi*x/2) + 0.5*pi*(-pi*sin(pi*x)*sin(pi*x/2) + 0.5*pi*cos(pi*x)*cos(pi*x/2)); % f(x)
h0 = 1; % h(0)
%g0 = -pi/2; % g(0)
g0 = 0; % dirichlet
x0 = 0; % starting point
err_stag = []; % norm err for standard grid
err_stan = []; % norm err for staggered grid
Ns = 2.^(5:12); % num grid points

for N = Ns
    dx = 1/N; % step size in space
    %x_stan = 0:dx:1; % includes boundary values, N+1 points
    %x_stan = dx:dx:1; % includes Neumann boundary value, N points
    x_stan = dx:dx:(1-dx); % dirichlet on both ends
    x_stag = dx/2:dx:(1-dx/2); % doesn't include boundary values, N-1 points
    
    [A1,f1] = discretize(0,h0,g0,N,dx,fx,d,x0); % standard grid
    u_stan = -A1\f1;
    err_stan = [err_stan; norm(u_stan'-usol(x_stan,0))];
    
    [A2,f2] = discretize(1,h0,g0,N,dx,fx,d,x0); % staggered grid
    [R,flag] = chol(-A2);
    if ~flag
        u_stag = R\(R'\f2);
    else
        disp('Factorization failed.')
        u_stag = -A2\f2;
    end
    err_stag = [err_stag; norm(u_stag'-usol(x_stag,0))];
    
    if N == 32
        
        figure(1); clf;
        plot(x_stan,u_stan,'.-'); hold on
        figure(2); clf;
        plot(x_stag,u_stag,'.-'); hold on
        
    elseif N == 64
        
        figure(1)
        plot(x_stan,u_stan,'.-')
        plot(x_stan,usol(x_stan,0)','-')
        xlabel('x')
        ylabel('u')
        title('Solution of steady state elliptic PDE on standard grid')
        legend('N=32','N=64','actual solution')
        
        figure(2)
        plot(x_stag,u_stag,'.-')
        plot(x_stag,usol(x_stag,0)','-')
        xlabel('x')
        ylabel('u')
        title('Solution of steady state elliptic PDE on staggered grid')
        legend('N=32','N=64','actual solution')
    end
    
end

figure(3); clf;
loglog(Ns,err_stan,'o-'); hold on
loglog(Ns,err_stag,'o-') % possibly change to dx instead of N
xlabel('N')
ylabel('||u(x,0)-u_{sol}(x,0)||')
title('Finite difference discretization')
legend('standard grid','staggered grid')
% % 1.3 hw5 non-linear heat equation
% using fd and (explicit) forward euler

d = @(x) 2+cos(2*pi*x); % d(x)
fxt = @(x,t) 0; % f(x,t)
ht = @(t) 0; % h(t)
gt = @(t) 0; % g(t), periodic domain
x0 = -1; % starting point
T = 10; % final time period
mu = 1; % CFL
N = 16; % grid resolution
dx = 2/N; % step size in space
dt = mu*(dx^2)/d(0); % step size in time

x_stan = (x0+dx):dx:(1-dx); % dirichlet at both ends, N-1 points

ux0 = @(u0,x) u0*(cos(0.5*pi*x)).^100;
u0 = 1;
un = ux0(u0,x_stan)';

fx = @(x)fxt(x,0);

[A,f] = discretize(2,ht(0),gt(0),N,dx,fx,d,x0); % standard grid

% set(gcf,'doublebuffer','on')
% v = VideoWriter('stag_cn.avi');
% open(v);

figure(1)
plot(x_stan,un)
drawnow()
%writeVideo(v,getframe)
pause(0.1)

for t = dt:dt:T

    %un1 = un + dt*A*un + dt*un.^2 + dt*f;
    un1 = un + (A+sparse(diag(un)))*dt*un + dt*f;
    un = un1;
    
    f(2:(N-2)) = fxt(x0+(2:(N-2))'*dx,t);
    f(N-1) = fxt(1-dx,t)+d(1-dx/2)*gt(t)/(dx^2);
    f(1) = fxt(x0+dx,t)+ht(t)*d(x0+dx/2)/(dx^2);

    figure(1)
    plot(x_stan',un)
    axis([-1 1 0 5])
    drawnow()
    %writeVideo(v,getframe)
    pause(0.1)

end
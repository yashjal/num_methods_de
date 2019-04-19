% 1.2 hw5 robustness of fd mol methods 
% for solving parabolic diffusion eq.

d = @(x) 1; % d(x)
fxt = @(x,t) 0; % f(x,t)
ht = @(t) 0; % h(t)
gt = @(t) 1; % g(t), assuming dirichlet
x0 = 0; % starting point
T = 0.01; % final time period
mu = 1; % CFL
N = 16; % grid resolution
dx = 1/N; % step size in space
dt = mu*(dx^2)/d(0); % step size in time

x_stan = dx:dx:(1-dx); % dirichlet at both ends, N-1 points
x_stag = dx/2:dx:(1-dx/2); % doesn't include boundary values, N-1 points

ux0 = @(x) double(x > 0.5); % discont initial cond
un_be = ux0(x_stan)';
un_cn = un_be;

fx = @(x)fxt(x,0);

[Ast,fst] = discretize(2,ht(0),gt(0),N,dx,fx,d,x0); % standard grid
fst1 = fst;

% set(gcf,'doublebuffer','on')
% v = VideoWriter('stag_cn.avi');
% open(v);

figure(1)
plot(x_stan,un_be); hold on
plot(x_stan,un_cn); hold off
drawnow()
%writeVideo(v,getframe)
pause(0.1)

for t = dt:dt:T

    fst1(2:(N-2)) = fxt(x0+(2:(N-2))'*dx,t);
    fst1(N-1) = fxt(1-dx,t)+d(1-dx/2)*gt(t)/dx^2;
    fst1(1) = fxt(x0+dx,t)+ht(t)*d(x0+dx/2)/(dx^2);

    un_be = be(Ast,un_be,fst1,dt);
    un_cn = cn(Ast,un_cn,fst,fst1,dt);

    fst = fst1;

    figure(1)
    plot(x_stan',un_be); hold on
    plot(x_stan',un_cn); hold off
    axis([0 1 0 1])
    drawnow()
    %writeVideo(v,getframe)
    pause(0.1)

end

sol = pdepe(0,@pdepde,@pdeic,@pdebc,(0:dx:1),(0:dt:T));
usol = sol(:,:,1);
figure(2)
plot((0:dx:1),usol(end,:)); hold on
plot(x_stan,un_be)
plot(x_stan,un_cn)
title('Testing robustness')
xlabel('x')
ylabel('u(x,T)')
legend('MATLAB soln','BE','CN')

function [c,f,s] = pdepde(x,t,u,dudx)
    c = 1;
    f = dudx;
    s = 0;
end

function u0 = pdeic(x)
    u0 = double(x > 0.5);
end

function [pl,ql,pr,qr] = pdebc(xl,ul,xr,ur,t)
    pl = ul;
    ql = 0;
    pr = ur-1;
    qr = 0;
end

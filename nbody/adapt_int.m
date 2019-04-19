G = 1; % gravitation constant
r = 1; % r(0) starting position of planet
dec = [1,0.5,0.25,0.125]; % decaying velocities for elliptical orbits

for d=dec
    
    v = d*r; % initial v(0), reduce for ellipse
    n = 2; % dimension of r
    N = 1; % number of planets
    m0 = 1; % mass of star
    r0 = zeros(n,1); % position of star

    r_0 = [r;0]; % intial position of planet
    v_0 = [0;v]; % inital velocity of planet
    E = -G*m0/r + 0.5*v^2; % energy
    a = -G*m0/(2*E); % semi-major axis
    T = 2*pi*a*(a/(G*m0))^0.5; % length of period
    tspan = [0 T];
    abs_tol = 10^(-3)*r; % absolute error tolerance
    y0 = [r_0; v_0]; % intital condition for modified first order system

    [t_sol,y_sol,delta_t,tall,xall] = bsrk23_ad(@(t,y) nbodyodef(t,y,n,G,m0,r0),tspan,y0,abs_tol);
    r_sol = y_sol(:,1:n*N); % computed solution
    rsol = r_0'; % actual solution at T

    disp(norm(r_sol(end,:)-rsol)); % norm error should meet target
    
    if d == 1
        figure(1);
        plot(tall,xall,'o');
        xlabel('t'); ylabel('x(t)');
        title('Adaptive integration for v(0)=1');
    else
        figure(2);
        plot(tall,xall,'o'); hold on;       
    end
    
    figure(3);
    plot(t_sol, delta_t, 'o-'); hold on;
    
end

figure(2);
xlabel('t'); ylabel('x(t)');
title('Adaptive integration');
legend('v(0)=0.5','v(0)=0.25','v(0)=0.125');

figure(3);
xlabel('t'); ylabel('\delta t');
title('Adaptive time steps');
legend('v(0)=0.5','v(0)=0.25','v(0)=0.125');
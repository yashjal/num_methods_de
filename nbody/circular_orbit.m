G = 1; % gravitation constant
r = 1; % length of planet position
w = 1; % angular frequency
n = 2; % dimension of r
N = 1; % number of planets
M = 22; % number of time steps
r0 = zeros(n,1); % position of star

T = 2*pi/w; % period of circular orbit 
m0 = r^3*w^2/G; % mass of star
delta_t = T/M;
rt = @(r,w,t)r*[cos(w*t);sin(w*t)]; % actual solution
r_0 = rt(r,w,0); % intial position of singular planet
v_0 = [0;w*r_0(1)]; % initial velcity of sing. planet

tspan = (0:delta_t:T)'; % interval of integration
rel_tol = 10^(-2)*r; % relative error tolerance
y0 = [r_0; v_0]; % intital condition for modified first order system

[t_sol,y_sol] = bsrk23(@(t,y) nbodyodef(t,y,n,G,m0,r0), tspan, y0);
r_sol = y_sol(:,1:n*N); % computed solution
rsol = zeros(length(tspan),n*N); % actual solution
for i=1:length(tspan)
   rsol(i,:) = rt(r,w,tspan(i))';
end

disp(abs(r_sol(end,1)-rsol(end,1))); %/abs(rsol(end,1)) for rel error
disp(abs(r_sol(end,2)-rsol(end,2))); %/abs(rsol(end,2)) for rel error
%[t_sol2,y_sol2] = ode23(@(t,y) nbodyodef(t,y,n,G,m0,r0), tspan, y0);

M2 = 2*M; % repeat for delta t/2
delta_t2 = T/M2;
tspan2 = (0:delta_t2:T)';

[t_sol2,y_sol2] = bsrk23(@(t,y) nbodyodef(t,y,n,G,m0,r0), tspan2, y0);
r_sol2 = y_sol2(:,1:n*N); % computed solution at time T

disp(abs(r_sol2(end,1)-rsol(end,1))); 
disp(abs(r_sol2(end,2)-rsol(end,2))); % abs error

r_sol3 = zeros(length(tspan),N*n); % solution for richardson extrapolation
p = 3; % bsrk23 is 3rd order accurate
cnt = 1;
%ref: http://bolyai.cs.elte.hu/~faragois/papers/stability_richardson_backward_euler_journal_version.pdf
for i = 1:2:length(tspan2)
    %r_sol3(cnt,:) = r_sol(cnt,:) + ((r_sol(cnt,:) - r_sol2(i,:))./(2^(p) - 1)) * delta_t^(p+1);
    r_sol3(cnt,:) = (2^(p)*r_sol2(i,:)-r_sol(cnt,:))/(2^(p)-1);
    cnt = cnt + 1;
end

disp(abs(r_sol3(end,1)-rsol(end,1)));
disp(abs(r_sol3(end,2)-rsol(end,2)));

% plot error vs time
figure(1); clf;
plot(tspan,vecnorm(r_sol-rsol,2,2),'ro-'); hold on;
plot(tspan,vecnorm(r_sol2(1:2:end,:)-rsol,2,2),'bo-'); hold on;
plot(tspan,vecnorm(r_sol3-rsol,2,2),'ko-'); hold on;
xlabel('t'); ylabel('e^{k}=||R^{K}-r(k\delta t)||_{2}');
title('Norm error in BS-RK3 solutions as time varies');
legend('\delta t=0.2856','\delta t/2','RE');

figure(2); clf;
semilogy(tspan,vecnorm(r_sol-rsol,2,2),'ro-'); hold on;
semilogy(tspan,vecnorm(r_sol2(1:2:end,:)-rsol,2,2),'bo-'); hold on;
semilogy(tspan,vecnorm(r_sol3-rsol,2,2),'ko-'); hold on;
xlabel('t'); ylabel('e^{k}=||R^{K}-r(k\delta t)||_{2}');
title('Norm error in BS-RK3 solutions as time varies');
legend('\delta t=0.2856','\delta t/2','RE');
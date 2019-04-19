G = 1; % gravitation constant
r = 1; % length of planet position
w = 1; % angular frequency
n = 2; % dimension of r
N = 1; % number of planets
M = 22; % number of time steps
r0 = zeros(n,1); % position of star

T = 2*pi/w; % period(s) of circular orbit 
m0 = r^3*w^2/G; % mass of star
delta_t = T/M;
T = 200 * T; % increase period number
rt = @(r,w,t)r*[cos(w*t);sin(w*t)]; % actual solution

tspan = (0:delta_t:T)'; % interval of integration
rel_tol = 10^(-2)*r; % relative error tolerance

rsol = zeros(length(tspan),n*N); % actual solution values
r_sol = zeros(length(tspan),n*N); % computed solution using VI
[t_sol,y_sol] = bsrk23(@(t,y) nbodyodef(t,y,n,G,m0,r0), tspan, y0);
r_sol2 = y_sol(:,1:n*N); % using BS-RK23

for i=1:length(tspan)
   rsol(i,:) = rt(r,w,tspan(i))';
   if i <= 2
       r_sol(i,:) = rsol(i,:);
   else
       r_sol(i,:) = F(1,r_sol(i-1,:)',n,G,m0,r0)' * delta_t^2 + 2*r_sol(i-1,:) - r_sol(i-2,:); % Verlet integrator
   end
end

disp(abs(r_sol(end,1)-rsol(end,1)));
disp(abs(r_sol(end,2)-rsol(end,2)));

figure(1); clf;
plot(tspan,vecnorm(r_sol-rsol,2,2),'ro-'); hold on;
plot(tspan,vecnorm(r_sol2-rsol,2,2),'go-');
xlabel('t'); ylabel('e^{k}=||R^{K}-r(k\delta t)||_{2}');
title('Norm error using Verlet Integrator and BS-RK23 over long period of time');
legend('Verlet Integrator','BS-RK23');

figure(2); clf;
plot(tspan,vecnorm(r_sol,2,2),'ro-'); hold on;
plot(tspan,vecnorm(r_sol2,2,2),'bo-');
plot(tspan,vecnorm(rsol,2,2),'ko-');
xlabel('t'); ylabel('||r(t)||');
title('Comparing solutions of BS-RK23 and Verlet');
legend('Verlet Integrator','BS-RK23','actual solution');
G = 1; % gravitation constant
n = 2; % dimension of planet position
N = 5; % number of planets
r0 = zeros(n,1); % position of star

m0 = 2; % mass of star
me = m0/332946; % mass of earth
m = [m0; 100*me; 2.5*me; 317.8*me; 1000*me; 10000*me]; % mass of planets
x_0 = [1; 1.25; 1.5; 1.8; 2.2]; % starting positions
r_0 = zeros(n*numel(x_0),1);
r_0(1:n:end) = x_0;
v_0 = zeros(n*numel(x_0),1); % initial velocities of planets
v_0(2:n:end) = 0.5*x_0; % assumes circular orbit

E = -G*m0/r_0(1) + 0.5*v_0(2)^2; % energy planet 1
a = -G*m0/(2*E); % semi-major axis planet 1
T = 2*pi*a*(a/(G*m0))^0.5; % T for planet 1
delta_t = 2.060769409233160e-04; % second min adaptive time stepping
T = 10*T;

tspan = (0:delta_t:T)'; % interval of integration
%tspan = [0 T];
y0 = [r_0; v_0];

[t_sol,y_sol] = bsrk23(@(t,y) nbodyodef(t,y,n,G,m,r0), tspan, y0);
%[t_sol2,y_sol2] = ode23(@(t,y) nbodyodef(t,y,n,G,m,r0), tspan, y0);

r_sol = y_sol(:,1:n*N);
xsol = r_sol(:,1:n:end);
ysol = r_sol(:,2:n:end);

%disp(diff(t_sol2)); % find min adaptive time steps

figure(1); clf;      % movie time!
loops = size(r_sol,1);
F(loops) = struct('cdata',[],'colormap',[]);
v = VideoWriter('yj627_solar_system.avi');
open(v);
for j = 1:100:loops
    scatter(0,0); hold on;
    scatter(xsol(j,:),ysol(j,:)); hold off;
    axis([-5 5 -5 5])
    drawnow
    F(j) = getframe;
    writeVideo(v,F(j));
end
close(v);
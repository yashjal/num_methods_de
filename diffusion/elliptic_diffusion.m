% Elliptic PDE of steady state diffusion eq. in 1D
% Solving for u(x,0) and estimating order of accuracy

d = @(x) 2 + cos(pi*x); % d(x)
usol = @(x,t) exp(-pi^2*t/4)*cos(0.5*pi*x); % manufactured soln
fx = @(x) 0.5*pi^2*cos(pi*x/2) + 0.5*pi*(-pi*sin(pi*x)*sin(pi*x/2) + 0.5*pi*cos(pi*x)*cos(pi*x/2)); % f(x)
h0 = 1; % h(0)
g0 = -pi/2; % g(0)
G = 0; % grid type, staggered or standard
N = 10000; % grid size
dx = 1/N; % step size in space
x_stan = 0:dx:1; % includes boundary values, N+1 points
x_stag = dx/2:dx:(1-dx/2); % doesn't include boundary values, N-1 points

if G == 0
    
    A = zeros(N+1,N+1); % non-symmetric because of Neumann BC
    f = zeros(N+1,1);
    for i = 2:N
        A(i,i-1) = d((i-1)*dx-dx/2);
        A(i,i+1) = d((i-1)*dx+dx/2);
        A(i,i) = -A(i,i-1)-A(i,i+1);
        f(i) = fx((i-1)*dx);
    end
    A(N+1,N) = A(N,N+1)+d(1+dx/2);
    A(N+1,N+1) = -A(N+1,N);
    A = A/dx^2;
    A(1,1) = 1;
    A = sparse(A);
    f(1) = -h0;
    f(N+1) = fx(1)+2*d(1+dx/2)*g0/dx;
    
    % solution u(x,0)
    u_G1 = -A\f;

end

if G == 1
    
    A = zeros(N,N); % snd in staggered grid
    f = zeros(N,1);
    for i = 2:N-1
        A(i,i-1) = d((i-1)*dx);
        A(i,i+1) = d(i*dx);
        A(i,i) = -A(i,i-1)-A(i,i+1);
        f(i) = fx((i-1)*dx+dx/2);
    end
    A(1,1) = -2*d(0)-A(2,1);
    A(1,2) = A(2,1);
    A(N,N-1) = A(N-1,N);
    A(N,N) = -A(N,N-1);
    A = A/dx^2;
    A = sparse(A);
    f(1) = fx(dx/2)+2*d(0)*h0/(dx^2);
    f(N) = fx(1-dx/2)+d(N*dx)*g0/dx;
    
    % solution u(x,0)
    [R,flag] = chol(-A);
    if ~flag
        disp('Factorization successful.');
        u_G2 = R\(R'\f);
    else
        disp('Factorization failed.');
        u_G2 = -A\f;
    end
end
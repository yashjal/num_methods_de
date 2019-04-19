function u = elliptic_solve(G,ht,gt,N,dx,fx,d)

% Elliptic PDE of steady state diffusion eq. in 1D
% Solving for u(x,0)

if G == 0 % standard drid
    
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
    f(1) = -ht;
    f(N+1) = fx(1)+2*d(1+dx/2)*gt/dx;
    
    % solution u(x,0)
    u = -A\f;

elseif G == 1 % staggered grid
    
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
    f(1) = fx(dx/2)+2*d(0)*ht/(dx^2);
    f(N) = fx(1-dx/2)+d(N*dx)*gt/dx;
    
    % solution u(x,0)
    [R,flag] = chol(-A);
    if ~flag
        disp('Factorization successful.');
        u = R\(R'\f);
    else
        disp('Factorization failed.');
        u = -A\f;
    end
    
end

end


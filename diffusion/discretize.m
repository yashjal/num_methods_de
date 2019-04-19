function [A,f] = discretize(G,ht,gt,N,dx,fx,d,x0)
% discretizes the 1D diffusion eq. in space

if G == 0 % standard grid
    
    A = zeros(N,N); % non-symmetric because of Neumann BC
    f = zeros(N,1);
    for i = 1:(N-1)
        if i ~= 1
            A(i,i-1) = d(x0+i*dx-dx/2);
            A(i,i+1) = d(x0+i*dx+dx/2);
            A(i,i) = -A(i,i-1)-A(i,i+1);
        else
            A(i,i+1) = d(x0+i*dx+dx/2);
            A(i,i) = -d(x0+dx/2)-A(i,i+1);
        end
        f(i) = fx(x0+i*dx);
    end
    A(N,N-1) = A(N-1,N)+d(1+dx/2);
    A(N,N) = -A(N,N-1);
    A = A/dx^2;
    A = sparse(A);
    f(1) = f(1)+ht*d(x0+dx/2)/dx^2;
    f(N) = fx(1)+2*d(1+dx/2)*gt/dx;

%     A = zeros(N+1,N+1); % non-symmetric because of Neumann BC
%     f = zeros(N+1,1);
%     for i = 2:N
%         A(i,i-1) = d(x0+(i-1)*dx-dx/2);
%         A(i,i+1) = d(x0+(i-1)*dx+dx/2);
%         A(i,i) = -A(i,i-1)-A(i,i+1);
%         f(i) = fx(x0+(i-1)*dx);
%     end
%     A(N+1,N) = A(N,N+1)+d(1+dx/2);
%     A(N+1,N+1) = -A(N+1,N);
%     A = A/dx^2;
%     A(1,1) = 1;
%     A = sparse(A);
%     f(1) = -ht;
%     f(N+1) = fx(1)+2*d(1+dx/2)*gt/dx;

elseif G == 1 % staggered grid
    A = zeros(N,N); % snd in staggered grid
    f = zeros(N,1);
    for i = 2:N-1
        A(i,i-1) = d(x0+(i-1)*dx);
        A(i,i+1) = d(x0+i*dx);
        A(i,i) = -A(i,i-1)-A(i,i+1);
        f(i) = fx(x0+(i-1)*dx+dx/2);
    end
    A(1,1) = -2*d(x0)-A(2,1);
    A(1,2) = A(2,1);
    A(N,N-1) = A(N-1,N);
    A(N,N) = -A(N,N-1);
    A = A/dx^2;
    A = sparse(A);
    f(1) = fx(x0+dx/2)+2*d(x0)*ht/(dx^2);
    f(N) = fx(1-dx/2)+d(N*dx)*gt/dx;
    
elseif G == 2 % standard grid with dirichlet on both boundaries
    
    A = zeros(N-1,N-1);
    f = zeros(N-1,1);
    for i = 1:(N-2)
        if i ~= 1
            A(i,i-1) = d(x0+i*dx-dx/2);
            A(i,i+1) = d(x0+i*dx+dx/2);
            A(i,i) = -A(i,i-1)-A(i,i+1);
        else
            A(i,i+1) = d(x0+i*dx+dx/2);
            A(i,i) = -d(x0+dx/2)-A(i,i+1);
        end
        f(i) = fx(x0+i*dx);
    end
    A(N-1,N-2) = A(N-2,N-1);
    A(N-1,N-1) = -A(N-1,N-2)-d(1-dx/2);
    A = A/dx^2;
    A = sparse(A);
    f(1) = f(1)+ht*d(x0+dx/2)/dx^2;
    f(N-1) = fx(1-dx)+d(1-dx/2)*gt/dx^2;
    
end

end


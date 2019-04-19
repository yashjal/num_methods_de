function u = cn(A,un,f1,f2,dt)
% crank-nicholson step in time
% f1 = fn, f2 = fn+1
I = speye(size(A));
X = (I-0.5*dt*A);
b = (I+0.5*dt*A)*un + dt*0.5*(f1+f2);

[R,flag] = chol(X);
if ~flag
    u = R\(R'\b);
else
    disp('Factorization failed.');
    u = -X\b;
end

end